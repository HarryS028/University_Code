# University coding projects

## Gauss-Seidel task 

Solve the elliptic PDE:

<img src="images%20for%20github/q2_equation.png" >

Subject to the boundary conditions:

<img src="images%20for%20github/q2_boundary.png" >

I first assumed:

<img src="images%20for%20github/q2_assumption.png" >

Where h is the grid spacing in the x-direction and k is the grid spacing in the y-direction. So our PDE becomes:

<img src="images%20for%20github/q2_grid.png" >

In our example, h and k are dx = 0.05 and dy = 0.05 respectively. As h=k we can rewrite the equations above as:

<img src="images%20for%20github/q2_final.png" >

Now we have an implicit scheme for solving elliptic PDE’s on a Cartesian grid, and can be solved by the Gauss_Seidel method in my code.

## Numerical_methods

This problem arose from the need to solve the following two equations modelling a two coupled pendula system. Here l is the length of 
each pendulum, k is the spring constant of the spring connecting the pendula and m1 and m2 are the masses of pendulum 1 and 2 
respectively. θ1 represents the angle of pendulum 1 to the vertical and θ2 represents the angle of pendulum 2 to the vertical. 

<img src="images%20for%20github/pendula_equations.png" >

Firstly, to linearise the set of equations I took a Taylor Series expansion of sin(θ) at 0, taking the linear terms I got that sin(θ) = θ as an approximation.  I then substituted θ into the original equations in place of sin(θ) to get the linearised equations. Seen below:

<img src="images%20for%20github/linearised_q1.png" >

As each equation is second order I split them up to create a system of first order ordinary differential equations to be solved simultaneously. 

Let θ1  =  θa  then:

<img src="images%20for%20github/q1_3.png" >

Also let θ2 = θx then:

<img src="images%20for%20github/q1_4.png" >

To apply the Adams Moulton method four initial values were needed for each of the variables I had used, in my c code I denoted θa, θb, θx, θy by thA, thB, thX, thY respectively. To generate these values I implemented a 4th order Runge-Kutta method at the start of my code. 

## Sputtering.mw 

This code uses Maple to model the affects and patterns created of ion sputtering on to a surface.
