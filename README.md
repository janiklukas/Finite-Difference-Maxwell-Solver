# Finite Difference Maxwell Solver

The functions in this project enable solving the Maxwell equations of electrodynamics (in potential form) for a number of point charges with fixed trajectories. The equation domain is a cube with periodic boundary conditions.

## Discretized Equations

We use the Lorenz gauge and Gaussian units with $c=1$. The governing equations for the potential fields $\phi$ and $\vec{A}$ are then given by:

$$(\partial_t^2-\nabla^2)\phi=4\pi\rho\quad\text{and}\quad(\partial_t^2-\nabla^2)\vec{A}=4\pi\vec{j}.$$

We divide the domain into $N^3$ equal cubic cells so that the Laplacian can be approximated as

$$(\nabla^2f)_ {i,j,k}\approx f_{i+1,j,k}+f_{i-1,j,k}+f_{i,j+1,k}+f_{i,j-1,k}+f_{i,j,k+1}+f_{i,j,k-1}-6f_{i,j,k}.$$

## Initial Field State

To compute the initial state we evaluate the static equations

$$ \nabla^2\phi=-4\pi\rho\quad\text{and}\quad\nabla^2\vec{A}=-4\pi\vec{j} $$

After discretization these can be recast into a $4N^3\times 4N^3$ linear system $Ax=b$ where $x$ contains the field values for all cells and $b$ the charge contributions. Since $A$ is a symmetric matrix, the system can be solved using the conjugate residual method with a sparse approximate inverse preconditioner.

## Time Evolution

A central difference method is used for the time derivatives:

$$\partial_t^2f(t)\approx\frac{1}{\Delta t^2}[f(t+\Delta t)-2f(t) + f(t-\Delta t)], $$

so two previous states are required at each time step.

This results in the following simple algorithm:

$$x^{(n+1)}=2x^{(n)}-x^{(n-1)}+\Delta t^2[Ax^{(n)}+b^{(n)}]$$

where $A$ is the same matrix as in the previous section.
