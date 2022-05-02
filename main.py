from src.assembly import *
from src.mesh import *
from src.solver import *
from src.visualization import *

def main():
    #testing
    f = lambda x: 1

    #creating mesh objects
    mesh = Mesh(4, delenuay=False)
    mesh_random = Mesh(6, randomize=True, sigma=0.2)
    mesh_refined = Mesh(20)

    #creating solvers.
    solver1 = Solver2D(f, mesh)
    solver2 = Solver2D(f, mesh_refined)
    solver3 = Solver2D(f, mesh_random)

    #solving the systems
    X1, Y1, U1 = solver1.solve()
    X2, Y2, U2 = solver2.solve()
    X3, Y3, U3 = solver3.solve()

    plot_mesh([mesh, mesh_random, mesh_refined], fname="figs/meshes.png",
        titles=[f"h = {round(mesh.h,3)}", f"h = {round(mesh_random.h,3)} at boundary", f"h = {round(mesh_refined.h,3)}"])

    plot_solution([[X1,Y1,U1],[X3,Y3,U3], [X2,Y2,U2]], fname="figs/solutions.png",
            titles=[f"h = {round(mesh.h,3)}", f"h = {round(mesh_random.h,3)} at boundary", f"h = {round(mesh_refined.h,3)}"])



if __name__ == "__main__":
    main()
