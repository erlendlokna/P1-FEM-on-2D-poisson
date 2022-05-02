from matplotlib import cm
import matplotlib.pyplot as plt

def plot_mesh(meshes, titles=None, size=(20, 5), fname=None):
    """
    meshes: array containing mesh objects
    titles: array containing titles (strings)
    fname: filename for saving
    """
    fig = plt.figure(figsize=size)

    for i in range(len(meshes)):

        ax = fig.add_subplot(1, len(meshes), i+1)

        mesh = meshes[i]

        if(titles): ax.set_title(titles[i])
        
        ax.triplot(mesh.X[:, 0], mesh.X[:, 1], mesh.T, c="gray")
        ax.scatter(mesh.X[:, 0], mesh.X[:, 1], s=20, c="black")

        ax.set_xlabel("x"); ax.set_ylabel("y")

        ax.grid(False)

    if(fname): plt.savefig(fname)
    plt.show()


def plot_solution(data, size=(18, 10), titles=None, cmap=True, fname=None):
    """
    data: array containing X, Y and U. Solutions to the equation.
    titles: array containing strings representing titles in the plot.
    """

    fig = plt.figure(figsize=size)

    for i in range(len(data)):

        ax = fig.add_subplot(1, len(data), i+1, projection='3d')

        X, Y, Z = data[i]

        if(cmap): ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.gnuplot, linewidth=0, antialiased=False)
        else: ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=False)

        if(titles): ax.set_title(titles[i])
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('u')

    if(fname): plt.savefig(fname)
    plt.show()