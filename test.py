from polyhedron_func import *
import numpy as np
from graph import *
import pdb
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import drawing

#
# V = np.array([[0,3],[3,0],[0,0],[3,3]])
# poly1 = polyhedron([V,np.array([])],'V')
# V = np.array([[1,1],[2,2],[2,1],[1,2]])
# poly2 = polyhedron([V,np.array([])],'V')
# poly12 = poly_intersect(poly1,poly2)
# poly12.minRep()
# # x = poly_push_out(poly,np.array([1.1,0.5]),0.5)
# H_set,sub_poly1,sub_poly2 = find_intersecting_hyperplane(poly1,poly2)
# pdb.set_trace()

g = Graph()
#
g.add_node("1")
g.add_node("2")
g.add_node("3")
g.add_node("4")
g.add_node("5")
g.add_node("6")

g.add_edge("1","2",1)
g.add_edge("1","4",2)
g.add_edge("1","3",1)
g.add_edge("2","3",2)
g.add_edge("3","4",1)
g.add_edge("3","5",2)
g.add_edge("4","6",1)
g.add_edge("5","6",2)
g.add_edge("2","5",1)

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


# g.remove_node("1")
# ss = g.reachable_node("2")
#
# V= np.array([[2.91182409, 3.2609008 ],
#        [2.9118006 , 3.2610172 ],
#        [3.96872476, 3.12465336],
#        [2.81472887, 3.76504068]])

# A = np.array([[ 0.04873852, -0.99881157],
#        [-0.12785426, -0.99179297],
#        [-0.98023608, -0.1978313 ],
#        [-0.4916411 ,  0.87079792],
#        [-0.98195434, -0.18911817]])
# b = np.array([-2.92751016, -3.60642758, -3.49938328,  1.89475321, -3.47597284])
#
#
# poly1 = polyhedron([V,np.array([])],'V')
# poly2 = polyhedron([A,b],'H')
# V,R = H2V(poly2.A,poly2.b)
# pos = np.array([4.02921651, 3.79685365])
# pdb.set_trace()
# poly=[]
# for i in range(0,4):
#     poly.append(polyhedron([V,np.array([])],'V'))
# for p in poly:
#     p.volume = 0.1
# pdb.set_trace()


# props = dict(boxstyle='round', facecolor='None',edgecolor='None')
# fig = plt.figure()
# ax = fig.gca()
#
#
#
# # ax.add_patch(map_patches[1])
# plt.axis('equal')
# plt.axis('off')
# plt.tight_layout()
#
# ax.add_patch(drawing.draw_convhull(poly1.V, ax, edgecolor='k', facecolor='y', alpha=0.2)[0])
# ax.add_patch(drawing.draw_convhull(poly2.V, ax, edgecolor='k', facecolor='y', alpha=0.2)[0])
#
#     # V=np.array([[0,0],[2,2],[2,0],[0,2]])
#     # ax.add_patch(pch.Polygon(V,closed=False,color='y',alpha=0.3))
# plt.show()
# H_set,sub_poly1,sub_poly2 = find_intersecting_hyperplane(poly1,poly2)
# pdb.set_trace()
