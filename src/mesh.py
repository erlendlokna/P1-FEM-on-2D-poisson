
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
import numpy as np

class Mesh:
    def __init__(self, axis_nodes, randomize=False, sigma=None, xlim=[0,1], ylim=[0,1], delenuay=True):
        """
        axis_nodes: nodes on each axis
        randomized: normal distribute inner nodes
        sigma: sigma value for normal distrubution
        delenuay: (bool) Delenuay triangulation
        
        """
        self.xlim = xlim #mesh limits.
        self.ylim = ylim
        self.Mx, self.My = axis_nodes, axis_nodes #number of nodes on axis
        self.M = axis_nodes * axis_nodes #number of nodes
        self.N = 2*(self.Mx-1)**2 #number of triangles.
        self.h = 1/(self.Mx-1)

        self.xs = np.linspace(self.xlim[0], self.xlim[1], self.Mx)
        self.ys = np.linspace(self.ylim[0], self.ylim[1], self.My)
        

        #getting X, 1d array of node coordinates.
        self.X = np.array([[i, j] for i in self.xs for j in self.ys]) #containg edge points/nodes.

        if(randomize): #using a normal distrubution on each coordinate.
            if(not sigma):
                sigma = 0.3 * self.h 

            for k in range(len(self.X)):
                i = k % self.Mx; j = k // self.My #the boundaries will not be randomized
                if(i==0 or i==self.Mx-1): continue
                if(j==0 or j == self.My-1): continue

                x = self.xlim[0]; y = ylim[0]

                while(True):
                    x = self.X[k][0] + self.X[k][0]*np.random.uniform(0, sigma)
                    y = self.X[k][1] + self.X[k][1]*np.random.uniform(0, sigma)
                    if((x > xlim[0] and x < xlim[1]) and (y > ylim[0] and y < ylim[1])): break

                self.X[k] = [x, y]

        #getting T, connectivity matrix
        if(delenuay or randomize): self.delaunayTriangulation()
        else: self.simpleTriangulation()

    def delaunayTriangulation(self):
        self.T = Delaunay(self.X).simplices
    
    def simpleTriangulation(self):
        #triangulates as done in class and in the notes of curry.
        self.T = []

        grid_indexes = np.reshape(np.arange(0, self.M, 1), (self.Mx, self.My))

        for i in range(len(grid_indexes)-1):
            current_row = grid_indexes[i]
            next_row = grid_indexes[i+1]

            for j in range(len(current_row)-1):
                self.T.append([current_row[j], current_row[j+1], next_row[j]])
                self.T.append([current_row[j+1], next_row[j+1], next_row[j]])