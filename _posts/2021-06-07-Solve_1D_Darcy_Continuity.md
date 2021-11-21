---
layout: post
title: Basic Python with Darcy's Law - Solving Darcy's Law and Continuity Equation on 1D Mesh
description: Create methods to solve for the pressure and velocity in a 1D saturated porous medium using Darcy's Law and the continuity equation
summary: Solving Darcy\'s law & the continuity equation over a mesh using different methods.
mathjax: true
---

While the code we\'ve developed so far allows us to go from parameters to plots in a few lines (assuming that all the classes and their methods are coded), the code also makes the assumption that the pressure is linear, which may not be true. To add more functionality to the code, let\'s figure out how to solve for the pressure distribution and velocity given only the parameters of the porous medium and the inlet and outlet conditions. 

### Math

If the pressure is also unknown, then we need another equation to solve for pressure. So far, we\'ve only used Darcy\'s Law. However, the assumptions of Darcy\'s law (slow flow through a porous medium) allow us to use another equation that accounts for mass conservation: the incompressible continuity equation:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\nabla\cdot\vec{u}&space;=&space;0"/>
</p>

Since our case is 1D, we can simplify this to:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\dfrac{\partial u}{\partial x}&space;=&space;0"/>
</p>

This means that the velocity in the pipe is constant in the x-direction, which we knew before. What we want to find is the pressure distribution so we need to plug in Darcy\'s law into the simplified equation above:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\dfrac{\partial u}{\partial x}&space;=&space;\dfrac{\partial}{\partial x}(-\dfrac{K}{\mu}\dfrac{dP}{dx})&space;=&space;0"/>
</p>

Our method from last week does not assume that the porous medium/liquid properties are constant so we\'ll keep the permeability and viscosity inside the parentheses. However, we can guess what the pressure distribution would be if we did. For now, let\'s move on to the equations that will be solved in Python.

In the previous [post](https://tim-munuhe.github.io/2021/05/14/mesh-fluidpm-plot), the velocity at the faces between the pressure points was coded but only briefly discussed. Let\'s flesh out the idea a bit more fully:

The mesh consists of individual cells that meet at their left and right faces, as in the figure below:

![Mid-mesh figure](../../../assets/images/Mid_mesh.png "Cells and faces in middle of mesh")

The middle cell experiences flow inwards from the left and outwards to the right. The adjacent cells have their own properties (K and mu specifically) so the pressure gradients within the cells may be different even though they have the same flow velocity (Darcy\'s Law). However, you also know that the pressure at the face is equal for both cell, so at face _i+0.5_:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?u_{i+0.5}&space;=&space;\dfrac{-K_i}{\mu_i}\dfrac{P_{i+0.5}-P_i}{\triangle x/2}&space;=&space;\dfrac{-K_{i+1}}{\mu_{i+1}}\dfrac{P_{i+1}-P_{i+0.5}}{\triangle x/2}"/>
</p>

We don\'t know the pressure at the face so we need to remove it so that we have an _equivalent_ permeability/viscosity factor that accounts for both cells. Luckily, we know K and mu in both cells so we can do some algebra to get the velocity across the face in terms of the adjacent cell pressures:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?u_{i+0.5}\;=\;-\dfrac{f_i%20f_{i+1}}{f_i+f_{i+1}}(P_{i+1}-P_i),\;\;\;\;\;f_i&space;=&space;%20\dfrac{K_{i}}{\mu_{i}}\dfrac{1}{x_{i+0.5}-x_i}"/>
</p>

This expression will take care of all interior faces. 

I wrote at the beginning that the continuity equation would allow use to solve for the pressure distribution, then combined it with Darcy\'s law. How does that work with our conception so far? Well, conservation of mass in a 1D mesh, assuming a constant cross-sectional area (i.e. all the faces have the same area), dictates that the velocity at both faces of cell _i_ be equal to each other:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?u_{i-0.5}\;=\;u_{i+0.5}"/>
</p>

We can combine this with the preceding equation to get our pressure relations. Since the cell of interest is cell _i_, we need _P<sub>i</sub>_ on one side and all the other terms on the other side:

 <p align="center">
<img src="https://latex.codecogs.com/svg.image?P_i\;=\;\dfrac{\dfrac{f_{i+1}f_i}{f_{i+1}+f_i}P_{i+1}+\dfrac{f_{i-1}f_i}{f_{i-1}+f_i}P_{i-1}} {\dfrac{f_{i+1}f_i}{f_{i+1}+f_i}+\dfrac{f_{i-1}f_i}{f_{i-1}+f_i}}"/>
</p>

At the boundaries, the pressures are given so the boundary cells have slightly different terms: the inner faces are interior faces so their factors don\'t change. For the outer (boundary) faces, the pressure at the face is known so the velocity calculation is easy there. For the most left cell (at origin), _P<sub>i</sub>_ is calculated as:

 <p align="center">
<img src="https://latex.codecogs.com/svg.image?P_i\;=\;\dfrac{\dfrac{f_{i+1}f_i}{f_{i+1}+f_i}P_{i+1}+f_iP_0}{\dfrac{f_{i+1}f_i}{f_{i+1}+f_i}+f_i}"/>
</p>

### Basic Solution using Gauss-Seidel

We have equations for the pressure of each cell in terms of each other, thereby creating a system of algebraic equations. If we had 3 or 4, we could solve them by hand by plugging in pressures we know for equations we don\'t. But that gets less efficient the finer the mesh becomes. That being said, the system of equations we have now allows us to debug in a straightforward manner by checking each individual calculation. Let\'s use that to our advantage now and use some more sophisticated methods later.

The most straightforward method I can think of to solve a system of algebraic equations is the [Gauss-Seidel method](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method). Skipping to the end, you substitute values of variables you do know forward into equations you do not. So, in our case, we know the inlet pressure (_P<sub>0</sub>_) so we can use that to find the pressure at the adjacent cell center _P<sub>1</sub>_. Notice though, that the pressure at that cell is also a function of the pressure at the next cell _P<sub>2</sub>_ that we do not know. So, we have a to make an initial guess then calculate this once or even multiple times. More reading is available at Gauss-Seidel's [wiki](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method) but some important points are:

* Gauss-Seidel is an iterative method requiring repeated solution of the same equations with changing variables
* Solutions (should[^GS]) improve with successive iterations until error is below some tolerance (convergence)
* The number of iterations it takes to convergence depends on the quality of the initial guess

So we\'ve prepped the math, let\'s jump to the code. First, the pseudocode:

```
# Solver Parameters
tol = 1E-6 # tolerance to determining stopping point of scheme
res = 1.0 # residual (initially greater than the tolerance
max_iter = 1000 # max iterations (so it doesn't go forever)
k = 0 # iteration counter


## Initialize
p[2:N-1] = zeros(N-2,1) # initial guess for cell centers

## Iteration Loop
do while ((res>tol)&&(k<max_iter))
  p_prev = p # previous iteration
  for i in range(1,N-1)
    Aw = ... # "west" factor (i-1)
	Ae = ... # "east" factor (i+1)
	Ap = Aw + Ae # center factor (i)
    p[i] = (Aw*p[i-1] + Ae*p[i+1])/Ap
  end
  res = sum(abs(p-p_prev)) # L2 norm of p_diff
  k += 1 # increase iteration count
end do
```

We first define the solver parameters: tol and max_iter define how soon the iteration loop terminates, depending either on convergence (tol) or a maximum iteration count. res and k represents how the scheme is doing and whichever reaches its limit sooner terminates the loop. Then, we initialize the middle mesh points' solution. Since we have pressure boundaries, p[0] and p[N-1] are already solved. Since we do not know how many iterations it will take to get a converged pressure, we use a while loop as the outer loop with the residual (defined here as the L2-norm of the change in the pressure vector per iteration) as the main escape condition. The iteration count is insurance so the program doesn't hang here.

Inside the loop, the Gauss-Seidel method allows the use of the next iteration's solution (p[i+1], p[i-1]) so the equations are are continually changing based on new data. There is another simple iterative method called the Jacobi method where the previous iteration's solution (p_prev in the pseudocode) is used to calculate the next iteration. The difference then becomes convergence speed, where Gauss-Seidel is generally faster. 

So what does this look like in Python? Not too different from the pseudocode but it is adjusted for better debugging:

```
def gauss_seidel(self,msh,pm): # need the mesh info and porous medium permeability

        # Solver Parameters
        tol = 1E-6 # tolerance to determining stopping point of scheme
        res = np.array([1.0],dtype=float) # residual (initially greater than the tolerance
        max_iter = 100 # max iterations (so it doesn't go forever)
        k = 0 # iteration counter

        # self.p[2:N-1] = zeros(N-2,1) # initial guess for cell centers
        p_samp = np.zeros([1,4],dtype=float)
        p_samp[0][:] = np.copy([self.p[1],self.p[3],self.p[msh.Nx-4],self.p[msh.Nx-2]])
        # p_samp[0][0] = np.copy(self.p[1])
        # p_samp[0][1] = np.copy(self.p[3])
        # p_samp[0][2] = np.copy(self.p[msh.Nx-4])
        # p_samp[0][3] = np.copy(self.p[msh.Nx-2])
        
        ## Iteration Loop
        while ((res[k]>tol)and(k<max_iter)):
            p_prev = np.copy(self.p)# previous iteration (copy to avoid using same mem loc)
            i = 1 # first cell center
            fw = -pm.K[i]/self.mu[i]/(msh.x[i]-msh.xc[i]) # f_i-1/2 -> f_i
            fe = -pm.K[i]/self.mu[i]/(msh.xc[i+1]-msh.x[i]) # f_i -> f_i+1/2
            fee = -pm.K[i+1]/self.mu[i+1]/(msh.x[i+1]-msh.xc[i+1]) #f_i+1/2 -> f_i+1
            Aw = fw
            Ae = fe*fee/(fe+fee)
            self.p[i] = (Aw*self.p[i-1] + Ae*self.p[i+1])/(Aw + Ae)
            for i in range(2,msh.Nx-2):
                fww = -pm.K[i-1]/self.mu[i-1]/(msh.xc[i]-msh.x[i-1]) # f_i-1->i-1/2
                fw = -pm.K[i]/self.mu[i]/(msh.x[i]-msh.xc[i]) # f_i-1/2 -> f_i
                fe = -pm.K[i]/self.mu[i]/(msh.xc[i+1]-msh.x[i]) # f_i -> f_i+1/2
                fee = -pm.K[i+1]/self.mu[i+1]/(msh.x[i+1]-msh.xc[i+1]) # f_i+1/2 -> f_i+1
                Aw = fw*fww/(fw+fww) # "west" factor (i-1 -> i)
                Ae = fe*fee/(fe+fee) # "east" factor (i -> i+1)
                self.p[i] = (Aw*self.p[i-1] + Ae*self.p[i+1])/(Aw + Ae)
            i = msh.Nx-2 # last cell center
            fww = -pm.K[i-1]/self.mu[i-1]/(msh.xc[i]-msh.x[i-1]) # f_i-1->i-1/2
            fw = -pm.K[i]/self.mu[i]/(msh.x[i]-msh.xc[i]) # f_i-1/2 -> f_i
            fe = -pm.K[i]/self.mu[i]/(msh.xc[i+1]-msh.x[i]) # f_i -> f_i+1/2
            Aw = fw*fww/(fw+fww) 
            Ae = fe 
            self.p[i] = (Aw*self.p[i-1] + Ae*self.p[i+1])/(Aw + Ae)
            p_samp = np.append(p_samp,[[self.p[1],self.p[3],self.p[msh.Nx-4],self.p[msh.Nx-2]]],axis=0)
            res = np.append(res,[sum(abs(self.p-p_prev))]) # L2 norm of p_diff
            k += 1 # increase iteration count
            # print(k,res[k-1])
            
        # Iterations are complete. Now for output
        print('Gauss-Seidel Complete. Iteration, Residual:',k,res[k])
        
        # I suggest using the pandas library for output to file. Compare the code below to the output
        # function coded from scratch in the mesh class
        res_vec = res[:,np.newaxis]
        df = pd.DataFrame(np.append(res_vec,p_samp,axis=1),columns=['res','x1','x3','x_N-4','x_N-2'])
        df.to_csv('GS_Out.csv',sep='\t')       
        return [k,res[k]]
```

The solver parameters remain the same but I\'ve added a pressure sampling array so that I can look at some exact numbers at different iterations. The coefficients Ae and Aw have also been constructed from factors representing what happens in halves of cells: fww represents the right half of the left cell, fw represents the left half of the center cell, and so on. The outermost cells have their pressures calculated outside of the loop, acknowledging that their outer faces are actually boundaries.  After the iteration loop, I report some data to the terminal (through print) and to a csv file through a pandas function. I hope you notice the much shorter code for outputting to a file with pandas compared to the scratch function in the mesh class. My recommendation: if you can\'t help fighting with your code, find tools that help you avoid those debugging battles. 

Since we\'re object-oriented, we can call this method quite easily using:
```
## Pressure calculation using Gauss-Seidel
fl_gs = fluid(base_mesh,base.fl)
print('Original P:',fl_gs.p[0:4],
      fl_gs.p[base_mesh.Nx-5:base_mesh.Nx-1])
>>> Original P: [0. 0. 0. 0.] [0. 0. 0. 0.]
[itera, res] = fl_gs.gauss_seidel(base_mesh,pm1)
print('Gauss-Seidel P:',fl_gs.p[0:4],
      fl_gs.p[base_mesh.Nx-5:base_mesh.Nx-1])
>>> Gauss-Seidel Complete. Iteration, Residual: 100 5.65338582708514
>>> Gauss-Seidel P: [ 0.         -0.00138227 -0.00462155 -0.00872932] [-80.27557742 -85.87430752 -91.51599336 -97.17199779]
```

Remembering our desired solution, it seems the scheme is way off. What\'s  going on? Remember, Gauss-Seidel is an iterative method, meaning it gets closer and closer to its solution with each iteration. We can run it 2 more times and compare the pressures to see this process:

```
[itera, res] = fl_gs.gauss_seidel(base_mesh,pm1)
print('Gauss-Seidel P:',fl_gs.p[0:4],
      fl_gs.p[base_mesh.Nx-5:base_mesh.Nx-1])
fl_gs2 = copy.deepcopy(fl_gs) # another 100 iterations
[itera, res] = fl_gs2.gauss_seidel(base_mesh,pm1)
fl_gs3 = copy.deepcopy(fl_gs2) # another 100 iterations
[itera, res] = fl_gs3.gauss_seidel(base_mesh,pm1)
>>> Plotting Code Here. Refer to Gauss_Seidel_Solve.py in Github Repo <<<
```

![Comparison Plot of P: True solution vs. GS.](../../../assets/images/p_comp.png "Pressure vs. x.")

The initial condition is the zero condition, which is good for the inlet but  far from the outlet solution. Since the solutions of individual nodes are determined by the neighbor, getting the final solution using this iterative method requires "communication" between the nodes, and it takes a while for this to occur.
So, at the very least, it is moving towards the correct solution, iteration by iteration. But this isn\'t an efficient way of solving it. Giving a better initial condition would definitely help and using better solution algorithms/methods can help. The successive over-relaxation (SOR) method takes a small step from the Gauss-Seidel method but, with some matrix set-up, I recommend using a direct solver.

In any case, our code/method works. It took some work to get here so let\'s try a problem that we don\'t know the answer to. Assume the viscosity changes linearly according to:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\mu\;=\;gx+b"/>
</p>

where _g_ and _b_ are constants. In truth, we can do some calculus and differential equation work to get the solution but, at first glance, it\'s not apparent what the solution should be. Even if you don\'t want to do the math, you can check the pressure solution by simply calculating the velocities at the faces to ensure they are constant.

```
class fluid():
.
.
.
    def mu_lin(self,mesh):
        N = mesh.Nx
        L = mesh.x[N-1]
        L0 = mesh.x[0]
        for i in range(1,N):
            self.mu[i] = (0.005-0.001)/(L-L0)*mesh.x[i]+0.001
.
.
.
fl_gs_mulin = fluid(base_mesh,base.fl)
print('Original P:',fl_gs_mulin.p[0:4],
      fl_gs_mulin.p[base_mesh.Nx-5:base_mesh.Nx-1])
fl_gs_mulin.p_lin(base_mesh)
print('Linear P:',fl_gs_mulin.p[0:4],
      fl_gs_mulin.p[base_mesh.Nx-5:base_mesh.Nx-1])
fl_gs_mulin.mu_lin(base_mesh)
for k in range(0,51): [itera, res] = fl_gs_mulin.gauss_seidel(base_mesh,pm1)
print('Final P:',fl_gs_mulin.p[0:4],
      fl_gs_mulin.mu[base_mesh.Nx-5:base_mesh.Nx-1])
>>> Original P: [0. 0. 0. 0.] [0. 0. 0. 0.]
>>> Linear P: [ 0. -1. -3. -5.] [-93. -95. -97. -99.]
.
.
.
>>> Gauss-Seidel Complete. Iteration, Residual: 1 9.391120880941628e-08
>>> Final P: [ 0.         -0.33986983 -1.05987099 -1.83320556] [0.00472 0.0048  0.00488 0.00496]
.
>>> Code to Plot mu, p, and u, also in the Gauss_Seidel_Solve.py file in Github Repo <<<
```

![mu and converged P for GS and linear viscosity.](../../../assets/images/node_conv_mu_P.png "mu/P vs. x.")

The pressure distribution is different now, but of course, the changing viscosity means that the pressure distribution must be higher in some places to maintain the same velocity throughout (remember, continuity equation/conservation of mass). This is reflected in the constant velocity profile shown below:

![converged u for GS and linear viscosity.](../../../assets/images/face_conv_u.png "u vs. x.")

In the next blog post, we\'ll use more efficient extant tools to solve this problem. 

#### Footnotes

[^GS]: To solve a linear system using the Gauss-Seidel method (or any of the methods I\'ll be discussing), the coefficient matrix needs to be diagonally dominant or symmetric positive definite. Otherwise, the code will not converge to a solution. This is a bit technical so just know, for the systems of equations I set up (where solutions at indivdual points are only calculated based on their neighbors), the matrices are always diagonally dominant. I urge caution before trying to use this on a linear system with random values or some "easy" linear system coming from a textbook. 