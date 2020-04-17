from scipy import sparse
import numpy as np
import cdd
from numpy.linalg import norm
from scipy.spatial import ConvexHull
import osqp
import pdb
from itertools import combinations


def H2V(A,b):
    # H representation to V representation
    dim = A.shape[1]
    n = A.shape[0]
    H=[]
    for i in range(0,n):
        Hi= [b[i]] + [-A[i][j] for j in range(0,dim)]
        H.append(Hi)
    H = cdd.Matrix(H)
    H.rep_type = cdd.RepType.INEQUALITY
    poly = cdd.Polyhedron(H)
    Vm = poly.get_generators()
    V=[]
    R=[]
    for i in range(0,Vm.row_size):
        if Vm[i][0]==1:
            V.append([Vm[i][j+1] for j in range(0,dim)])
        else:
            R.append([Vm[i][j+1] for j in range(0,dim)])
    return np.array(V),np.array(R)

def V2H(V,R):
    # V representation to H representation
    dim = V.shape[1]
    nv = V.shape[0]
    nr = R.shape[0]
    Vm=[]
    sv = np.ones([nv,1])

    V = np.concatenate((sv,V),axis=1)
    if nr>0:
        sr = np.zeros([nr,1])
        R = np.concatenate((sr,R),axis=1)
        VR = np.concatenate((V,R))
    else:
        VR = V
    VR = cdd.Matrix(VR)
    VR.rep_type = cdd.RepType.GENERATOR
    poly = cdd.Polyhedron(VR)
    Hm = poly.get_inequalities()

    A=[]
    b=[]
    for i in range(0,Hm.row_size):
        b.append(Hm[i][0])
        A.append([-Hm[i][j+1] for j in range(0,dim)])
    return np.array(A),np.array(b)
class polyhedron():
    def __init__(self,invar,type):
        if type =='H':
            self.A = invar[0]
            self.b = invar[1]
            self.V, self.R = H2V(self.A,self.b)
            self.dim = self.A.shape[1]
            self.volume = []
            self.isempty = False
            if self.V.size == 0:
                self.isempty = True
            if not self.isempty:
                self.center = np.mean(self.V,axis=0)
        elif type == 'V':
            self.V = invar[0]
            self.R = invar[1]
            self.A, self.b = V2H(self.V,self.R)
            self.dim = self.V.shape[1]
            self.volume = []
            self.isempty = False
            if self.V.size == 0:
                self.isempty = True
            if not self.isempty:
                self.center = np.mean(self.V,axis=0)
        else:
            error("Type not supported, must be either 'V' or 'H'")
    def contain(self,p):
        # return True or False, and the margin of violation
        if not isinstance(p,np.ndarray):
            p = np.array(p)
        return all(np.matmul(self.A,p)<=self.b),min(self.b-np.matmul(self.A,p))
    def poly_contain(self,poly):
        # if a polytope is contained in self
        if not isinstance(p,polyhedron):
            return False
        res = True
        for v in poly.V:
            if not self.contain(v)[0]:
                res = False
                break
        return res
    def get_volume(self):
        if not self.volume:
            if self.R:
                self.volume = np.inf
            else:
                volume = ConvexHull(self.V).volume
        return self.volume
    def get_facet_vertices(self):
        # return the index of vertices contained in each facet
        facets = [[] for i in range(0,self.A.shape[0])]
        for i in range(0,self.V.shape[0]):
            axb = self.b-np.matmul(self.A,self.V[i])
            for j in range(0,axb.shape[0]):
                if abs(axb[j])<1e-6:
                    # pdb.set_trace()
                    facets[j].append(self.V[i])

        return facets
    def minRep(self):
        # get the minimum representation
        n = self.V.shape[0]
        V_rmv_idx = []
        for i in range(0,n):
            if min(self.b-np.matmul(self.A,self.V[i]))>1e-5:
                V_rmv_idx.append(i)
        self.V = np.delete(self.V,V_rmv_idx,axis=0)
        if self.V.shape[0]==0:
            self.isempty = True
        else:
            A_rmv_idx = []

            for i in range(0,self.A.shape[0]):
                try:
                    Avb = np.tile(self.b[i],(1,n))-np.matmul(self.A[i],self.V.transpose())
                except:
                    pdb.set_trace()
                if min(Avb[0])>1e-5:
                    A_rmv_idx.append(i)
            self.A = np.delete(self.A,A_rmv_idx,axis=0)
            self.b = np.delete(self.b,A_rmv_idx)

def poly_intersect(poly1,poly2):
    if poly1.dim!=poly2.dim:
        error("dimension do not match!")
    else:
        A = np.concatenate((poly1.A,poly2.A)).copy()
        b = np.concatenate((poly1.b,poly2.b)).copy()
        poly12 = polyhedron([A,b],'H')
        poly12.minRep()
        return poly12


def poly_projection(poly,x,margin = 0):
    # project a point x to a polyhedron with margin, margin>0 means the projection is on the shrinked polyhedron
    dim = x.shape[0]
    # pdb.set_trace()
    if poly.contain(x)[0]:
        return x,0.0

    P = sparse.csc_matrix(np.eye(dim))
    q = -x

    # pdb.set_trace()
    A = sparse.csc_matrix(poly.A)
    l = -np.inf*np.ones(poly.b.shape)
    u = np.array(poly.b)
    for i in range(0,poly.b.shape[0]):
        u[i] = u[i]-margin*norm(poly.A[i])

    prob = osqp.OSQP()
    prob.setup(P, q, A, l, u, alpha=1.0,verbose=False)
    res = prob.solve()
    x_proj = res.x
    if res.info.status_val == 1 or res.info.status_val == 2:
        dis = norm(x_proj-x)
        return x_proj,dis
    else:
        return x,0.0

def minkowski_shrink_poly(poly,radius):
    # shrink or bloat(radius<0) a polyhedron with minkowski sum
    A = np.array(poly.A)
    b = np.array(poly.b)

    for i in range(0,b.shape[0]):
        b[i] = b[i]-radius*norm(A[i])
    new_poly = polyhedron([A,b],'H')
    return new_poly
def Hausdorff_distance(poly1,poly2):
    # Hausdorff distance of two polyhedron (or polyunion, in which case it's pseudo Hausdorff distance)
    if poly1.isempty:
        return [],[]
    if isinstance(poly2,polyhedron):
        dis_list = []
        for i in range(0,poly1.V.shape[0]):
            x_proj,dis = poly_projection(poly2,poly1.V[i])
            dis_list.append(dis)
        return max(dis_list),min(dis_list)
    elif isinstance(poly2,list) and isinstance(poly2[0],polyhedron):
        dis_list = []
        for i in range(0,poly1.V.shape[0]):
            dis_V = []
            for poly in poly2:
                x_proj,dis = poly_projection(poly,poly1.V[i])
                dis_V.append(dis)
            dis_list.append(min(dis_V))
        return max(dis_list),min(dis_list)
    elif isinstance(poly2,dict) and isinstance(list(poly2.values())[0],polyhedron):
        dis_list = []
        for i in range(0,poly1.V.shape[0]):
            dis_V = []
            for poly in poly2.values():
                x_proj,dis = poly_projection(poly,poly1.V[i])
                dis_V.append(dis)
            dis_list.append(min(dis_V))
        return max(dis_list),min(dis_list)
def dis2bdry(poly,x):
# distance of a point to the boundary of poly, if x is insde, distance >0, otherwise <0
    if isinstance(poly,polyhedron):
        dis_list = []
        axb = poly.b-np.matmul(poly.A,x)
        for i in range(0,poly.A.shape[0]):
            dis_list.append(axb[i]/norm(poly.A[i]))
        return min(dis_list),[]
    elif isinstance(poly,list) and isinstance(poly[0],polyhedron):
        dis_list = []
        for poly_i in poly:
            dis_list_i = []
            axb = poly_i.b-np.matmul(poly_i.A,x)
            for i in range(0,poly_i.A.shape[0]):
                dis_list_i.append(axb[i]/norm(poly_i.A[i]))
            dis_list.append(min(dis_list_i))
        dis_list = np.array(dis_list)
        return max(dis_list),dis_list.argmax()
    elif isinstance(poly,dict) and isinstance(list(poly.values())[0],polyhedron):
        dis_list = []
        for idx in poly.keys():
            dis_list_i = []
            axb = poly[idx].b-np.matmul(poly[idx].A,x)
            for i in range(0,poly[idx].A.shape[0]):
                dis_list_i.append(axb[i]/norm(poly[idx].A[i]))
            dis_list.append(min(dis_list_i))
        dis_list = np.array(dis_list)
        return max(dis_list),list(poly.keys())[dis_list.argmax()]
def poly_push_out(poly,x,radius):
# find a point that is at least radius away from poly that is closed to x, nonconvex, so no guarnatee on global optimality
    A_norm = np.zeros(poly.A.shape[0])
    for i in range(0,poly.A.shape[0]):
        A_norm[i]=norm(poly.A[i])
    axb = poly.b-np.matmul(poly.A,x)+radius*A_norm
    if min(axb)<0:
        return x
    else:
        idx = axb.argmin()
        x = x+poly.A[idx]*axb[idx]/(A_norm[idx]**2)
        return x
def find_intersecting_hyperplane(poly1,poly2):
    # For two intersecting polyhedrons, get subsets that are irreducable on both sides, i.e., smallest possible subpolys
    # that satisfies subpoly1 \cup subpoly2 \cup (poly1\cap poly2)=poly1 \cup poly2
    poly12 = poly_intersect(poly1,poly2)
    if poly12.isempty is True:
        return [],[],[]
    bdry_pts = []
    non_bdry_pts = []

    n = poly12.V.shape[0]
    axb1 = np.tile(poly1.b,(n,1)).transpose()-np.matmul(poly1.A,poly12.V.transpose())
    axb2 = np.tile(poly2.b,(n,1)).transpose()-np.matmul(poly2.A,poly12.V.transpose())
    for i in range(0,n):
        if min(axb1[:,i])<=1e-9 and min(axb2[:,i])<=1e-9:
            bdry_pts.append(poly12.V[i])
        else:
            non_bdry_pts.append(poly12.V[i])
    bdry_pts = np.array(bdry_pts)
    m = bdry_pts.shape[0]
    H_set = []
    sub_poly1 = []
    sub_poly2 = []
    bdry_idx = range(0,len(bdry_pts))
    if poly12.dim == 2:
        comb = combinations(bdry_idx, 2)
    elif poly12.dim == 3:
        comb = combinations(bdry_idx, 3)

    for idx in list(comb):
        if poly12.dim == 2:
            x1 = bdry_pts[idx[0]]
            x2 = bdry_pts[idx[1]]
            if norm(x1-x2)>1e-4:
                A_orth = (x2-x1)/norm(x1-x2)
                A = np.array([-A_orth[1],A_orth[0]])
                b = np.matmul(A,x1)
        elif poly12.dim == 3:
            x1 = bdry_pts[idx[0]]
            x2 = bdry_pts[idx[1]]
            x3 = bdry_pts[idx[2]]
            if norm(x1-x2)>1e-4 and norm(x1-x3)>1e-4 and norm(x3-x2)>1e-4:
                A = np.cross(x1-x2, x2-x3)
                A = A/norm(A)
                b = np.matmul(A,x1)

        exist = False
        for hyperplane in H_set:
            if norm(hyperplane[0]-A)<1e-5 and norm(hyperplane[1]-b)<1e-5:
                exist = True
                break
            elif norm(hyperplane[0]+A)<1e-5 and norm(hyperplane[1]+b)<1e-5:
                exist = True
                break
        if exist is False:
            Avb = (np.tile(b,(1,m))-np.matmul(A,bdry_pts.transpose()))[0]

            other_idx = list(range(0,len(bdry_pts)))
            for i in idx:
                other_idx.remove(i)

            if all(Avb[other_idx]<=0) or  all(Avb[other_idx]>=0):

                if all(Avb>=0):
                    A = -A
                    b = -b

                H_set.append([A,b])
                n1 = poly1.V.shape[0]
                n2 = poly2.V.shape[0]
                belong = "poly2"
                for i in range(0,n1):
                    if b-np.matmul(A,poly1.V[i])>0:
                        if not poly2.contain(poly1.V[i])[0]:
                            belong = "poly1"
                            break
                if belong == "poly1":
                    AA = np.array([A])
                    A_bar = np.concatenate((poly1.A,AA))
                    b_bar = np.append(poly1.b,b)
                    sub_poly1.append(polyhedron([A_bar,b_bar],'H'))
                elif belong == "poly2":
                    AA = np.array([A])
                    A_bar = np.concatenate((poly2.A,AA))
                    b_bar = np.append(poly2.b,b)
                    sub_poly2.append(polyhedron([A_bar,b_bar],'H'))
    return H_set,sub_poly1,sub_poly2
