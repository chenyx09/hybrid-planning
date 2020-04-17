import numpy as np
import irispy
import osqp
import json
import cdd
import pdb
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as pch
from matplotlib import animation
from collections import defaultdict
from numpy.linalg import norm
from polyhedron_func import *
from graph import *
import drawing
import sys, getopt
import logging
from tulip import transys, spec, synth, dumpsmach




def ray_projection(x0,poly,ray,radius):
    # "fake" lidar function
    A1 = np.matmul(poly.A,ray)
    b1 = poly.b-np.matmul(poly.A,x0)
    p_idx = np.where(A1>1e-6)
    n_idx = np.where(A1<-1e-6)
    zero_idx = np.where(abs(A1)<1e-6)
    if min([1]+b1[zero_idx].tolist())<0:
        return radius

    else:
        t_lb = max([0]+(b1[n_idx]/A1[n_idx]).tolist())
        t_ub = min([radius]+(b1[p_idx]/A1[p_idx]).tolist())

        if t_lb>t_ub:
            return radius
        else:
            return t_lb

class map():
    def __init__(self,lb,ub,obs_or_infile,type='V'):
        self.lb = lb
        self.ub = ub
        self.dim = len(lb)
        self.obs_set = []
        if self.dim ==2:
            N = 64   #lines of lidar rays
            self.lidar_theta = np.linspace(0,2*np.pi/N*(N-1),N)
            self.lidar_ray = []
            for i in range(0,N):
                self.lidar_ray.append([np.cos(self.lidar_theta[i]),np.sin(self.lidar_theta[i])])
            self.lidar_ray = np.array(self.lidar_ray)
        if isinstance(obs_or_infile, str):
            with open(obs_or_infile, 'r') as filehandle:
                obs_or_infile = json.load(filehandle)

        for obs in obs_or_infile:
            self.obs_set.append(polyhedron([obs[0],obs[1]],type))
    def plot_map(self,ax):
        patches = []
        if self.dim==2:
            patches.append(pch.Rectangle(self.lb,self.ub[0]-self.lb[0],self.ub[1]-self.lb[1],\
            0.0,edgecolor='k',fill=False))
            for obs in self.obs_set:
                patches.extend(drawing.draw_convhull(obs.V, ax, edgecolor='k', facecolor='b', alpha=0.2))
            return patches
        elif self.dim==3:
            print("not written yet")
        else:
            error("Not visualizable")
    def lidar_scan(self,x0,lidar_radius,ap): # ap means atomic proposition, i.e., area of interest
        if not (all(x0>=self.lb) and all(x0<=self.ub)):
            return []

        point_cloud = []
        point_type = []
        ap_points = []
        ap_idx = []
        for i in range(0,len(self.lidar_ray)):

            min_d = lidar_radius
            if self.lidar_ray[i][0]>1e-6:
                min_d = min(min_d,(self.ub[0]-x0[0])/self.lidar_ray[i][0])
            elif self.lidar_ray[i][0]<-1e-6:
                min_d = min(min_d,(self.lb[0]-x0[0])/self.lidar_ray[i][0])

            if self.lidar_ray[i][1]>1e-6:
                min_d = min(min_d,(self.ub[1]-x0[1])/self.lidar_ray[i][1])
            elif self.lidar_ray[i][1]<-1e-6:
                min_d = min(min_d,(self.lb[1]-x0[1])/self.lidar_ray[i][1])
            lidar_point = x0+self.lidar_ray[i]*min_d+np.random.random(self.dim)*1e-3
            for obs in self.obs_set:
                t = ray_projection(x0,obs,self.lidar_ray[i],lidar_radius)
                if t<min_d:
                    min_d = t
                    lidar_point = x0+t*self.lidar_ray[i]+np.random.random(self.dim)*1e-3
            point_cloud.append(lidar_point)
            if min_d == lidar_radius:
                point_type.append(1)
            else:
                point_type.append(0)

            for p in ap.keys():
                t = ray_projection(x0,ap[p],self.lidar_ray[i],lidar_radius)
                if t<min_d:
                    lidar_point = x0+t*self.lidar_ray[i]
                    ap_points.append(lidar_point)
                    ap_idx.append(p)
        return point_cloud,point_type,ap_points,ap_idx
    def check_collision(self,pos,radius):
        collision = False
        obs_collision = []
        for obs in self.obs_set:
            x_proj,dis=poly_projection(obs,pos)
            if dis<radius:
                collision = True
                obs_collision = obs
                break
        return collision, obs_collision
class way_point():
    def __init__(self,pos,dis,free_space_idx):
        self.pos = pos
        self.dis = dis
        self.free_space_idx = free_space_idx
class robot_agent():
    def __init__(self,pos,env,radius=0.1,dt=0.1,v_max=1.0,lidar_radius=3.0):
        self.pos = pos                            # position
        self.env = env                            # map
        self.dim = pos.shape[0]
        self.radius = radius                      # radius of the robot
        self.dt = dt                              # sampling time
        self.v_max = v_max                        # maximum velocity
        self.free_space = {}                      # polytopes that represents the free space
        self.lidar_radius = lidar_radius          # range of lidar
        self.pot_wp = []                          # potential way points
        self.wp = []                              # way points
        self.discrete_graph = Graph()             # discrete graph in terms of the free space
        self.free_space_edge = defaultdict(list)  # position of the intersections of free spaces
        self.explored_wp = []                     # way points already visited
        self.t = 0
        self.path = []                            # discrete state path
        self.ap = {}                              # atomic propositions, including the name as keys and set as values
        self.label = defaultdict(list)            # label function: free space -> atomic proposition
        self.refine_success = True                # flag that the free space refinement is successful
        self.reachable_ap = []                    # atomic propositions whose set is reachable from curent position
        self.Tulip_spec = None                    # LTL spec for Tulip
        self.discrete_controller = None           # symbolic controller to be synthesized by Tulip
        self.prior_discrete_graph = []            # discrete graph when symbolic controller was last synthesized
        self.prior_label = []                     # label function when symbolic controller was last synthesized
        self.mode = 'explore'                     # mode, explore or patrol
        self.discrete_ctrl_cmd = []               # control command from the symbolic controller

    def step(self,vel):
        # single integrator dynamics, can be modified
        self.pos = self.pos + vel*self.dt
        self.t = self.t+self.dt
    def add_ap(self,name,poly):
        # add an atomic proposition with name and set as a polytope
        self.ap[name]=poly
    def add_Tulip_spec(self,spec):
        # add Tulip specification
        self.Tulip_spec = spec
    def vel_to_wp(self,waypoint):
        # calculate the desired velocity to the way point
        v_des = (waypoint-self.pos)/self.dt
        if norm(v_des)>self.v_max:
            return v_des/norm(v_des)*self.v_max
        else:
            return v_des
    def extend_free_space(self,point_cloud,point_type,ap_points,ap_idx):
        # use lidar scan to add new free spaces
        obstacles = grouping_point_cloud(point_cloud)
        potential_set = []
        if self.free_space:
            idx = self.determine_discrete_state()
            dis_set = []
            for pt in point_cloud:     # choose the point in the point cloud that is far away from existing free space
                dis_set.append(dis2bdry(self.free_space,pt)[0]+dis2bdry(self.free_space[idx],pt)[0])
            dis_set = np.array(dis_set)
            selected_pt = point_cloud[dis_set.argmin()]
            print(selected_pt)
        else:
            selected_pt = self.pos

        bounds = irispy.Polyhedron.from_bounds(self.env.lb, self.env.ub)
        region = irispy.inflate_region(obstacles, 0.5*(self.pos+selected_pt), bounds=bounds, return_debug_data=False)
        # call IRIS to compute the free space, with center chosen as the mid point between current position and the point selected from the point cloud
        new_poly = polyhedron([region.polyhedron.getA(),region.polyhedron.getB()],'H')

        # shrink the polytope with radius of the robot
        new_poly_shrinked = minkowski_shrink_poly(new_poly,1.1*self.radius)
        if new_poly_shrinked.contain(self.pos)[0]:
            potential_set.append(new_poly_shrinked)

        # repeat the same procedure with the center being the current position
        region = irispy.inflate_region(obstacles, self.pos, bounds=bounds, return_debug_data=False)
        new_poly = polyhedron([region.polyhedron.getA(),region.polyhedron.getB()],'H')

        new_poly_shrinked = minkowski_shrink_poly(new_poly,1.1*self.radius)
        potential_set.append(new_poly_shrinked)

        # if any ap currently not reachable can be seen by the lidar, extend free space towards that direction
        checked_ap = list(self.reachable_ap)
        for i in range(0,len(ap_points)):
            if ap_idx[i] not in checked_ap:
                region = irispy.inflate_region(obstacles, ap_points[i], bounds=bounds, return_debug_data=False)
                new_poly = polyhedron([region.polyhedron.getA(),region.polyhedron.getB()],'H')

                new_poly_shrinked = minkowski_shrink_poly(new_poly,1.1*self.radius)
                checked_ap.append(ap_idx[i])
                if new_poly_shrinked.contain(self.pos)[0]:
                    potential_set.append(new_poly_shrinked)

        # for all potential free spaces, remove the unnecessary intersections with existing free spaces and add some of them
        for set in potential_set:
            if not self.free_space:
                self.add_new_free_set(set)

            else:
                set_to_add = [set]
                for fs in self.free_space.values():
                    for new_set in set_to_add:

                        H_set,sub_poly1,sub_poly2 = find_intersecting_hyperplane(fs,new_set)
                        if sub_poly2:
                            set_to_add.remove(new_set)
                            set_to_add = set_to_add + sub_poly2

                for new_set in set_to_add:
                    added = False
                    for prop in self.ap.keys():
                        if prop not in self.reachable_ap:
                            poly_sect = poly_intersect(self.ap[prop],new_set)
                            if not poly_sect.isempty:  # if the new set intersects with an ap set that was not reachable, add it
                                if not added:
                                    # seperate out the intersection, and add the original new set
                                    self.add_new_free_set(new_set)

                                self.add_new_free_set(poly_sect)
                                idx = list(self.free_space.keys())[list(self.free_space.values()).index(poly_sect)]
                                self.label[idx].append(prop)
                                self.update_pot_wp(1.2*self.radius)

                    if Hausdorff_distance(new_set,self.free_space)[0]>3*self.radius:
                        # if there is substantial progress measured in Hausdorff distance, add it
                        if not added:
                            self.add_new_free_set(new_set)
                            self.update_pot_wp(1.2*self.radius)

        for pt in point_cloud:
            # add new potential waypoints
            dis,free_space_idx = dis2bdry(self.free_space,pt)
            if dis<=-self.radius:
                # the potential way points are the projections of point cloud onto the free space
                wp_projected = poly_projection(self.free_space[free_space_idx],pt,self.radius)[0]
                self.pot_wp.append(way_point(wp_projected,dis,free_space_idx))
        self.update_reachable_ap()

    def update_reachable_ap(self):
        # update the atomic propositions reachable from current position
        idx = self.determine_discrete_state()
        self.reachable_ap=[]
        reachable_set = self.discrete_graph.reachable_node(idx)
        for fs in reachable_set:
            for ap in self.label[fs]:
                if ap not in self.reachable_ap:
                    self.reachable_ap.append(ap)

    def determine_discrete_state(self):
        # determine the index of the free space of the current robot position
        if not self.free_space:
            return []
        else:
            candidate = []
            dis2cen = []
            for idx in self.free_space.keys():
                if self.free_space[idx].contain(self.pos)[1]>-1e-2:
                    candidate.append(idx)
                    dis2cen.append(norm(self.pos-self.free_space[idx].center))
            if not candidate:
                pdb.set_trace()
                return []

            else:
                dis2cen = np.array(dis2cen)
                return candidate[dis2cen.argmin()]
    def refine_free_space(self,point_cloud,point_type,idx):
        # when observing lidar points violating the free space, refine the free space
        added = False
        pts = [pt for pt,type in zip(point_cloud,point_type) if type==0]

        violated_pts = []
        # bloat the polytope back
        bloated_poly = minkowski_shrink_poly(self.free_space[idx],-1.1*self.radius)

        max_violation = 0
        for pt in pts:
            # pdb.set_trace()
            iscontain, dis = bloated_poly.contain(pt)
            if iscontain:
                violated_pts.append(pt)
                max_violation = max(max_violation,dis)
        if max_violation>0.1*self.radius:
            obstacles = grouping_point_cloud(violated_pts)
            facets = bloated_poly.get_facet_vertices()

            for facet in facets:
                if facet:
                    # generate obstacles based on the facets of the previous polytope
                    m = self.dim+1-len(facet)
                    m = max(1,m)
                    new_obs = np.concatenate((np.array(facet),np.array(facet[0]+1e-4*np.random.random((m,self.dim)))))
                    obstacles.append(new_obs.transpose())



            bounds = irispy.Polyhedron.from_bounds(self.env.lb, self.env.ub)
            region = irispy.inflate_region(obstacles, self.pos, bounds=bounds, return_debug_data=False)

            deleted_fs = self.free_space[idx]
            self.delete_free_space(idx)

            new_poly = polyhedron([region.polyhedron.getA(),region.polyhedron.getB()],'H')

            new_poly_shrinked = minkowski_shrink_poly(new_poly,1.1*self.radius)
            if new_poly_shrinked.R.shape[0]>0:
                # returns unbounded polytopes, need to debug
                pdb.set_trace()


            if not self.free_space:
                self.add_new_free_set(new_poly_shrinked)

            else:
                set_to_add = [new_poly_shrinked]
                for fs in self.free_space.values():
                    for new_set in set_to_add:

                        H_set,sub_poly1,sub_poly2 = find_intersecting_hyperplane(fs,new_set)
                        if sub_poly2:
                            set_to_add.remove(new_set)
                            set_to_add = set_to_add + sub_poly2


                for new_set in set_to_add:
                    if new_set.R.shape[0]>0:
                        pdb.set_trace()
                    if Hausdorff_distance(new_set,self.free_space)[0]>2*self.radius or new_set.contain(self.pos)[0]:
                        if not added:
                            added = True
                        self.add_new_free_set(new_set)
                        self.update_pot_wp(1.2*self.radius)
            if added:
                self.select_wp()
                self.update_reachable_ap()
                if not self.wp:
                    # no waypoint is returned, need debug
                    pdb.set_trace()
            else:
                # if no set is added, the refinement fails, then add back the original polytope and try next time step
                self.add_new_free_set(deleted_fs)
                self.update_pot_wp(1.2*self.radius)


            return added
        else:
            # no violation, return []
            return []

    def update_pot_wp(self,threshold):
        # update potential waypoints based on the free space, remove all points with distance to the free space below threshold
        for wp in self.pot_wp:
            dis,free_space_idx = dis2bdry(self.free_space,wp.pos)
            wp.free_space_idx = free_space_idx
            if dis >=threshold:
                self.pot_wp.remove(wp)
            elif dis<self.radius:
                wp.pos = poly_projection(self.free_space[free_space_idx],wp.pos,self.radius)[0]
    def sparse_wp(self,threshold=0.1):
        # remove visited waypoints that are too close to others
        selected_wp = []
        for pt in self.explored_wp:
            select = True
            for pt1 in selected_wp:
                if norm(pt-pt1)<threshold:
                    select = False
                    break
            if select:
                selected_wp.append(pt)
        self.explored_wp = selected_wp
    def select_wp(self):
        # select waypoints from potential ones based on scores
        self.wp = []
        self.path = []
        if not self.pot_wp:

            point_cloud,point_type,ap_points,ap_idx = self.env.lidar_scan(self.pos,self.lidar_radius,self.ap)
            self.extend_free_space(point_cloud,point_type,ap_points,ap_idx)
            if not self.pot_wp:
                print("No waypoint of interest found")
                return []
        start_idx = self.determine_discrete_state()
        # select waypoint from reachable free spaces
        reachable_set = self.discrete_graph.reachable_node(start_idx)
        wp_score = []
        self.sparse_wp(threshold=self.radius)
        for wp in self.pot_wp:
            if wp.free_space_idx in reachable_set:
                min_dis = 100

                for wp1 in self.explored_wp:
                    min_dis = min(min_dis,norm(wp.pos-wp1))
                # score = -(minimum distance to visited waypoint) - (distance to the boundary(negative when outside the free space))- distance to current position
                wp_score.append(-wp.dis+2*min_dis-0.01*norm(wp.pos-self.pos))
            else:
                wp_score.append(-np.inf)
        wp_score = np.array(wp_score)
        selected_wp = self.pot_wp[wp_score.argmax()]
        if selected_wp.free_space_idx not in reachable_set:
            return False
        for wp in self.pot_wp:
            if norm(wp.pos-selected_wp.pos)<2*self.radius:
                self.pot_wp.remove(wp)
        dest_idx = selected_wp.free_space_idx

        if isinstance(start_idx,list):
            # error, need debug
            pdb.set_trace()
            # dijsktra search to get discrete path
        visited, path = dijsktra(self.discrete_graph, start_idx,dest_idx)
        if not path:
            return False
        for i in range(0,len(path)-1):
            # add the middle waypoints between the free spaces along the path
            self.wp.append(self.free_space_edge[(path[i],path[i+1])])

        wp_candidate,dis = poly_projection(self.free_space[dest_idx],selected_wp.pos,self.radius)

        print(wp_candidate)

        self.wp.append(wp_candidate)
        self.path = path
        return True
    def add_edge_free_space(self,idx1,idx2):
        # add an edge between free spaces with an nonempty intersection
        assert(idx1 in self.free_space.keys())
        assert(idx2 in self.free_space.keys())

        poly12 = poly_intersect(self.free_space[idx1],self.free_space[idx2])
        if poly12.isempty:
            if (idx1,idx2) in self.free_space_edge.keys():
                del self.free_space_edge[(idx1,idx2)]
            if (idx2,idx1) in self.free_space_edge.keys():
                del self.free_space_edge[(idx2,idx1)]
            self.discrete_graph.remove_edge(idx1,idx2)
            return False
        else:
            mid_point = np.mean(poly12.V,axis=0)
            self.free_space_edge[(idx1,idx2)] = mid_point
            self.free_space_edge[(idx2,idx1)] = mid_point
            dist12 = norm(self.free_space[idx1].center-mid_point)+norm(self.free_space[idx2].center-mid_point)
            self.discrete_graph.add_edge(idx1,idx2,dist12)
            return True

    def add_new_free_set(self,new_free_set):
        # add a new free space
        idx = 0
        while idx in self.free_space.keys():
            idx = idx+1

        self.free_space[idx]=new_free_set
        self.discrete_graph.add_node(idx)
        for i in self.free_space.keys():
            self.add_edge_free_space(i,idx)
    def delete_free_space(self,idx):
        # remove a free space (mostly due to the refinement)
        self.discrete_graph.remove_node(idx)
        # pdb.set_trace()
        del self.free_space[idx]
        # pdb.set_trace()
        for entry in list(self.free_space_edge.keys()):
            if entry[0] == idx or entry[1] == idx:
                del self.free_space_edge[entry]

    def tulip_synthesis(self):
        # call Tulip to synthesize a discrete controller
        sys = transys.FTS()
        sys.states.add_from(self.free_space.keys())
        idx = self.determine_discrete_state()
        sys.states.initial.add(idx)

        for node in self.discrete_graph.edges.keys():
            sys.transitions.add_comb(set([node]), set(self.discrete_graph.edges[node]))

        sys.atomic_propositions.add_from(set(self.ap.keys()))

        for state in self.label.keys():
            sys.states.add(state,ap=set(self.label[state]))

        ctrl = synth.synthesize(self.Tulip_spec, sys=sys)
        if ctrl is None:
            pdb.set_trace()
        else:
            dumpsmach.write_python_case("discrete_ctrl.py", ctrl, classname="symb_ctrl")
            from discrete_ctrl import symb_ctrl
            self.discrete_controller = symb_ctrl()
            self.prior_discrete_graph = self.discrete_graph.copy()
            self.prior_label = self.label.copy()



    def explore_controller_step(self):
        # controller when mode = explore
        collision, obs_collision = self.env.check_collision(self.pos,self.radius)
        if collision:
            wp_found = False
            if self.wp:
                wp_candidate = self.wp[0]
            else:
                wp_candidate = self.pos
            collision, obs_collision = self.env.check_collision(wp_candidate,self.radius)
            counter = 0
            while collision:

                wp_candidate = poly_push_out(obs_collision,wp_candidate,self.radius*1.2)
                collision, obs_collision = self.env.check_collision(wp_candidate,self.radius)
                counter = counter+1
                if counter>5:
                    print("cannot find collision free way point")
            self.wp = [wp_candidate]
            self.path = [None]
            pdb.set_trace()
            point_cloud = []
            point_type = []

        else:
            point_cloud,point_type,ap_points,ap_idx = self.env.lidar_scan(self.pos,self.lidar_radius,self.ap)
            idx = self.determine_discrete_state()
            if idx:
                dis = dis2bdry(self.free_space[idx],self.pos)[0]
                # pdb.set_trace()
                if (dis>self.radius*0.01 and int(self.t/self.dt)%5==0) or not self.refine_success:
                    self.refine_success = self.refine_free_space(point_cloud,point_type,idx)

            if int(self.t/self.dt)%5==0:
                self.explored_wp.append(self.pos)
            if not self.wp:
                self.extend_free_space(point_cloud,point_type,ap_points,ap_idx)

                success = self.select_wp()
                if not success:
                    pdb.set_trace()
                    wp = np.random.rand(self.dim)
                    fs_idx = self.determine_discrete_state()
                    wp,dis = poly_projection(self.free_space[fs_idx],wp)
                    self.wp = [wp]
                print(self.pos)

        vel = self.vel_to_wp(self.wp[0])

        self.step(vel)
        if norm(self.pos-self.wp[0])<self.radius*0.2:
            # reaching a waypoint
            self.explored_wp.append(self.wp[0])
            self.wp.remove(self.wp[0])
            self.path.remove(self.path[0])
        return point_cloud,point_type

    def LTL_controller_step(self):
        # controller when mode = patrol
        collision, obs_collision = self.env.check_collision(self.pos,self.radius)

        if collision:
            wp_found = False
            wp_candidate = self.pos
            counter = 0
            while collision:

                wp_candidate = poly_push_out(obs_collision,wp_candidate,self.radius*1.2)
                collision, obs_collision = self.env.check_collision(wp_candidate,self.radius)
                counter = counter+1
                if counter>5:
                    print("cannot find collision free way point")
            self.wp = [wp_candidate]
            self.path = [None]
            pdb.set_trace()
            point_cloud = []
            point_type = []

        else:
            point_cloud,point_type,ap_points,ap_idx = self.env.lidar_scan(self.pos,self.lidar_radius,self.ap)
            idx = self.determine_discrete_state()

            dis = dis2bdry(self.free_space[idx],self.pos)[0]
            if ((dis>self.radius*0.01 and int(self.t/self.dt)%5==0)) or self.refine_success is False:
                self.refine_success = self.refine_free_space(point_cloud,point_type,idx)
                if self.refine_success is True:
                    if not graph_similar(self.discrete_graph,self.prior_discrete_graph)\
                    or not self.label==self.prior_label:
                        self.mode = 'explore'

            if not self.wp:
                # Let the finite state controller take one step to get the next discrete state
                self.discrete_ctrl_cmd = self.discrete_controller.move()
                next_idx = self.discrete_ctrl_cmd['loc']
                while idx==next_idx:
                    self.discrete_ctrl_cmd = self.discrete_controller.move()
                    next_idx = self.discrete_ctrl_cmd['loc']

                if next_idx in self.discrete_graph.edges[idx]:
                    self.path = [next_idx]
                    self.wp.append(self.free_space_edge[(idx,next_idx)])
                else:
                    visited, path = dijsktra(self.discrete_graph, idx,next_idx)
                    if not path:
                        return False
                    for i in range(0,len(path)-1):
                        self.wp.append(self.free_space_edge[(path[i],path[i+1])])
                    self.path = path[0:-1]

                print(self.pos)

        vel = self.vel_to_wp(self.wp[0])
        self.step(vel)
        if norm(self.pos-self.wp[0])<self.radius*0.2:

            self.wp.remove(self.wp[0])
            self.path.remove(self.path[0])
        return point_cloud,point_type

def simulate_robot(robot,T,record=False,outfile=[]):
    dt = robot.dt
    nframe = int(T/dt)
    fig = plt.figure()
    ax = plt.gca()

    def animate(t,robot):
        print(t*robot.dt)
        if set(robot.ap.keys()).issubset(set(robot.reachable_ap)) and robot.mode is 'explore':
            robot.tulip_synthesis()
            robot.mode = 'patrol'

        if robot.mode is 'explore':
            point_cloud,point_type = robot.explore_controller_step()
        elif robot.mode is 'patrol':
            point_cloud,point_type = robot.LTL_controller_step()
        ax.clear()
        map_patches = robot.env.plot_map(ax)

        for patch in map_patches:
            ax.add_patch(patch)
        plt.axis('equal')
        plt.axis('off')
        plt.tight_layout()

        plt.plot(robot.pos[0],robot.pos[1],'mX',markersize=8)
        for i in range(0,len(point_cloud)):
            if point_type[i]==0:
                plt.plot(point_cloud[i][0],point_cloud[i][1],'go',markersize=3)
            elif point_type[i]==1:
                plt.plot(point_cloud[i][0],point_cloud[i][1],'ro',markersize=3)
        for wp in robot.wp:

            plt.plot(wp[0],wp[1],'ko')

        props = dict(boxstyle='round', facecolor='None',edgecolor='None')
        for idx in robot.free_space.keys():

            ax.text(robot.free_space[idx].center[0], robot.free_space[idx].center[1], str(idx),  fontsize=10,
        verticalalignment='top', bbox=props)
            ax.add_patch(drawing.draw_convhull(robot.free_space[idx].V, ax, edgecolor='k', facecolor='y', alpha=0.2)[0])

        for prop in robot.ap.keys():

            ax.text(robot.ap[prop].center[0], robot.ap[prop].center[1]+0.5, str(prop),  fontsize=10,
        verticalalignment='top', bbox=props)
            ax.add_patch(drawing.draw_convhull(robot.ap[prop].V, ax, edgecolor='k', facecolor='r', alpha=0.5)[0])

        ax.text(3, 5, str(robot.mode),  fontsize=13,
    verticalalignment='top', bbox=props)


    anim = animation.FuncAnimation(fig, animate, fargs=(robot,),
                                   frames=nframe,
                                   interval=dt*1000,
                                   blit=False, repeat=False)
    plt.show()
    if record:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        anim_name = outfile
        anim.save(anim_name,writer=writer)



def grouping_point_cloud(point_cloud,threshold = 0.2):
    dim = len(point_cloud[0])
    pts = []
    obstacles = []
    last_point = point_cloud[0]
    i = 0
    while i<len(point_cloud):
        if norm(point_cloud[i]-last_point)<threshold:
            last_point = point_cloud[i]
            pts = pts+[last_point]
            i = i+1
            if len(pts)>=dim+1:
                obstacles.append(np.array(pts).transpose())
                pts = []
                if i<len(point_cloud):
                    last_point = point_cloud[i]

        else:
            m = dim+1-len(pts)

            pts = np.concatenate((np.array(pts),np.random.random((m,dim))*0.2*threshold+last_point))
            obstacles.append(pts.transpose())
            i = i+1
            pts = []
            if i<len(point_cloud):
                last_point = point_cloud[i]
    if len(pts)>0:
        m = dim+1-len(pts)
        pts = np.concatenate((np.array(pts),np.random.random((m,dim))*0.2*threshold+pts[0]))
        obstacles.append(pts.transpose())
    return obstacles

def main(argv):
    record = False
    outfile = ''
    try:
        opts, args = getopt.getopt(argv,'hr:o:')
    except getopt.GetoptError:
        print('Robot_Exploration.py -r <True or False> -o <ofile>')
        pdb.set_trace()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Robot_Exploration.py -r <True or False> -o <ofile>')
            sys.exit()
        elif opt in ("-r", "--record"):
            if arg in ('True','1','true','t','yes'):
                record = True
            else:
                record = False
        elif opt in ("-o","--outfile"):
            outfile = arg
    if record:
        matplotlib.use("Agg")
        if not outfile:
            outfile = 'robot_animation.mp4'
    obs = [[[1,0],[6,0],[6,1],[1,1]],\
          [[0,2],[2,2],[2,3],[0,3]],\
          [[3,2],[5,2],[5,3],[3,3]],\
          [[0,4],[4,4],[4,5],[0,5]],\
          [[5,4],[6,4],[6,5],[5,5]]]
    for i in range(0,len(obs)):
        obs[i]=[np.array(obs[i]),np.array([])]
    lb=[0,0]
    ub=[6,5]
    env = map(lb,ub,obs,'V')
    V1 = np.array([[4.3,4.3],[4.3,4.7],[4.7,4.3],[4.7,4.7]])
    poly1 = polyhedron([V1,np.array([])],'V')
    V2 = np.array([[0.3,3.3],[0.3,3.7],[0.7,3.3],[0.7,3.7]])
    poly2 = polyhedron([V2,np.array([])],'V')

    x0 = np.array([0.5,0.5])


    robot = robot_agent(x0,env)
    ap = {'dest1':poly1,'dest2':poly2}
    env_vars = {}
    env_init = set()
    env_prog = set()
    env_safe = set()

    sys_vars = {}
    sys_init = set()
    sys_prog ={'dest1','dest2'}
    sys_safe = set()
    specs = spec.GRSpec(env_vars, sys_vars, env_init, sys_init,
                        env_safe, sys_safe, env_prog, sys_prog)

    specs.moore = True
    specs.qinit = '\E \A'
    robot.add_Tulip_spec(specs)
    for p in ap.keys():
        robot.add_ap(p,ap[p])


    simulate_robot(robot,80,record=record,outfile=outfile)


if __name__ == '__main__':
    main(sys.argv[1:])
