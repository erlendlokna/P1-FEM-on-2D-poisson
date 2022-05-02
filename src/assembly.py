import numpy as np
from src.mesh import *

def add_submatrix(A, A_new, t):
    #function for matrix operation. Returns "A[t,t] + A_new" where t is indices to verticies of the triangle beeing evaluated.
    A_out = A.copy()

    for i, t_i in enumerate(t):
        for j, t_j in enumerate(t):
            A_out[int(t_i)][int(t_j)] += A_new[j][i] #matrix operation.

    return A_out

def add_subarray(F, F_new, t):
    #function for operation "F[t] + F_new", where t is a list of verticies.
    F_out = F.copy()

    for i, t_i in enumerate(t):
        F_out[int(t_i)] += F_new[i]
    
    return F_out

def quadrature2D(f, t):
    """
    lam = np.array([
        [.333333333333333, .333333333333333, .333333333333333],
        [.081414823414554, .459292588292723, .459292588292723],
        [.658861384496480, .170569307751760, .170569307751760],
        [.898905543365938, .050547228317031, .050547228317031],
        [.008394777409958, .263112829634638, .728492392955404]
    ])
    
    rho = np.array([.144315607677787, .095091634267285, .103217370534718, .032458497623198, .027230314174435])
    """
    """
    lam=np.array([[.666666666666667, .166666666666667, .166666666666667]])
    rho=np.array([.333333333333333])
    """
    lam=np.array([[1/3, 1/3, 1/3]]) #using the simplest nodes and weights
    rho = np.array([1])

    F = np.zeros(3)

    for i in range(len(lam)):
        x_lam = lam[i, 0] * t[0] + lam[i, 1] * t[1] + lam[i, 2] * t[2]
        for j in range(3):
            F[j] = F[j] + rho[i] * lam[i, j] * f(x_lam)
    return F


class Assembly:
    def __init__(self, f, mesh):
        #assembles the mesh. creates stiffness matrix and F vector.

        self.mesh = mesh
        self.f = f

        self.A = np.zeros((mesh.M, mesh.M))
        self.F = np.zeros(mesh.M)

        rhs = np.array([[1, 0], [0, 1], [-1, -1]])

        for k in range(mesh.N): #looping through each triangle in the mesh.
            t = mesh.T[k] #triangle indexes to verticies.

            triangle = np.array([mesh.X[int(t[0])], mesh.X[int(t[1])], mesh.X[int(t[2])]]) #verticies
            
            J = np.array([triangle[0] - triangle[2], triangle[1] - triangle[2]])

            G = np.matmul(rhs, np.linalg.inv(J))

            abs_det_J = np.abs(np.linalg.det(J))
            
            A_new = 0.5 * np.matmul(G, np.transpose(G)) * abs_det_J
            F_new = 0.5 * abs_det_J * quadrature2D(self.f, triangle)

            self.A = add_submatrix(self.A, A_new, t)
            self.F = add_subarray(self.F, F_new, t)

        self.add_boundaries() #adding the boundary conditions to A and F.
        
    def add_boundaries(self):
        #adds dirichlet boundary conditions to the stiffness matrix and load vector.
        #At boundary points, x_i, this function adds identity rows to A[i][j] and 0 to F[i].
        #this ensures that u=0 at the boundaries.
        for i in range(self.mesh.Mx):
            self.A[i] = np.zeros(self.A.shape[0]); self.A[i][i] = 1
            self.F[i] = 0

        #top boundaries:
        for i in range(self.A.shape[0]-self.mesh.Mx, self.A.shape[0]):
            self.A[i] = np.zeros(self.A.shape[0]); self.A[i][i] = 1
            self.F[i] = 0

        #side boundaries:
        grid_indexes = np.reshape(np.arange(0, self.mesh.Mx * self.mesh.My, 1), (self.mesh.Mx, self.mesh.My))
        for i in grid_indexes[1:-1, 0]: #left
            self.A[i] = np.zeros(self.A.shape[0]); self.A[i][i] = 1
            self.F[i] = 0
        for i in grid_indexes[1:-1, -1]: #right
            self.A[i] = np.zeros(self.A.shape[0]); self.A[i][i] = 1
            self.F[i] = 0    