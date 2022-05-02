from src.assembly import *
import numdifftools as nd
from scipy.interpolate import interpolate

class Solver2D(Assembly):
    #solver for the 2D FEM poisson
    def __init__(self, f, mesh):
        super().__init__(f, mesh)        
    
    def solve(self):
        #solves the 2d poisson
        self.X, self.Y = np.meshgrid(self.mesh.xs, self.mesh.ys)
        self.Uh = np.reshape(np.linalg.solve(self.A, self.F), (self.mesh.Mx, self.mesh.My)) #numnpy meshgrid format.
        self.u_interp = interpolate.interp2d(self.X, self.Y, self.Uh, kind='linear') #interpolation function
        
        return self.X, self.Y, self.Uh
    
    def interp(self, x, y):
        #finds a interpolation estimate between the calculated node values
        return self.u_interp(x, y)[0]

    def errorL2(self, u):
        # Calculates the L2 error by integrating over each triangle.

        error2 = 0

        for k in range(len(self.mesh.T)):
            #looping through the triangles

            t = self.mesh.T[k] #triangle

            triangle = np.array([self.mesh.X[int(t[0])], self.mesh.X[int(t[1])], self.mesh.X[int(t[2])]]) #vertex coordinates

            def _f(x): return (self.interp(x[0], x[1]) - u(x[0], x[1]))**2

            error2 += np.sum(quadrature2D(_f, triangle))

        return np.sqrt(error2)
    
    def errorH1(self, u):
        #calculates the H1 norm estimate.
        
        error2 = 0

        for k in range(len(self.mesh.T)):
            #looping through the triangles

            t = self.mesh.T[k] #triangle

            triangle = np.array([self.mesh.X[int(t[0])], self.mesh.X[int(t[1])], self.mesh.X[int(t[2])]]) #vertex coordinates

            def _f(x): return (nd.Gradient(u)(x[0], x[1]) - nd.Gradient(self.interp)(x[0], x[1]))**2

            error2 += np.sum(quadrature2D(_f, triangle))

        return np.sqrt(self.errorL2(u) + error2)
