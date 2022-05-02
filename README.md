
This scipt solves numerically the 2D poisson equation presented below. It uses a finite element method on a square grid. The triangulation is done using a Delenuay algorithm. 

2D poisson equation:

<img src="https://render.githubusercontent.com/render/math?math=\nabla u = u_{xx} + u_{yy} = f(x,y)">

However, this is a quick implementation and creates a solution that resembles an exact solution. The error is not yet checked
