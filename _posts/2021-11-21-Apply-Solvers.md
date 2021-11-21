---
layout: post
title: Basic Python with Darcy's Law - Algorithms to Solve System of Linear Equations
description: Implement algorithms to solve system of linear equations created from Darcy\'s Law and Continuity Equation
summary: Solving Darcy\'s law & the continuity equation over a mesh using commonly used custom and packaged algorithms.
mathjax: true
---

We have equations for the pressure at every point within the 1D domain, but we're using a basic method to solve the equations. Here's some trivia: Gauss mentioned it first in 1823 but the method wasn't published until Seidel wrote it up in 1874. Since then, many newer methods have been developed to solve linear systems of equations. This post will cover three more commonly used algorithms: Cholesky decomposition, Gauss-Jordan, and the conjugate gradient method.

## System Characterization

First, let's consider our system of equations: to find the pressure, we solve only three unique sets of equations: the equations for the cells near the inlet and outlet and the equation for all the cells in the middle. Additionally, the cells' pressures only depend on the adjacent pressures, of which two are given (inlet and outlet). Such a system can be described in matrix form Ax = b, where A is a coefficient matrix, x is a column vector representing the variables we want to find, and b is a column vector representing the constants in each equation:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?A&space;=&space;\begin{bmatrix}1&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;0&space;\\A_{21}&space;&&space;-&space;A_{21}-&space;A_{23}&space;&&space;A_{23}&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;0&space;\\...&space;&&space;A_{32}&space;&&space;-A_{32}&space;-&space;A_{34}&space;&&space;A_{34}&space;&&space;...&space;&&space;...&space;&&space;0&space;\\...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;\\...&space;&&space;...&space;&&space;A_{i,i-1}&space;&&space;-A_{i,i-1}&space;-&space;A_{i,i&plus;1}&space;&&space;A_{i,i&plus;1}&space;&&space;...&space;&&space;...&space;\\...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;\\0&space;&&space;&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;...&space;&&space;1&space;\\\end{bmatrix}"/>
</p>

<p align="center">

<img src="https://latex.codecogs.com/svg.image?x&space;=&space;\begin{bmatrix}P_1&space;\\:&space;\\P_i&space;\\:&space;\\P_N\end{bmatrix},&space;\&space;b&space;=&space;\begin{bmatrix}P_{inlet}&space;\\:&space;\\0&space;\\:&space;\\P_{outlet}\end{bmatrix}&space;"/>

</p>

The inlet and outlet pressures are represented by the top and bottom rows of the matrix and both vectors, respectively. Because they have known pressures, they only have one value in the matrix and the given pressure in the b-vector. All other "internal" pressures depend on their neighboring pressures only. Therefore, their b-vector values are zero and they have three entries in their A-matrix (the left, center, and right cells). 

One important feature of the matrix: it is symmetric. That is, that the coefficient matrix's transpose is equal to the matrix itself. We can see this instantly by checking the entries on either side of the top left to bottom right diagonal. Additionally, the matrix is positive definite. This is not as obvious, unless you've had prior experience with such a matrix, but it can be checked by conducting Gaussian elimination to get the matrix into upper triangular/row echelon form. Let\'s get the row-echelon form using Python and NumPy:

```
import numpy as np

def row_reduc(A_org):
    N,M = A_org.shape
    if (N != M): raise Exception('Matrix must be square! Size:',N,',',M)
    
    # Setup
    A = np.copy(A_org)
    A_inv = np.eye(N) # A inverse initialized as identity matrix
    P_gen = np.eye(N) # General Permutation matrix
    P = np.eye(N) # Permutation matrix for this column pivot

    # Reduce to upper tridiagonal coefficient matrix
    for i in range(0, N-1):
        # Find row of pivot element in column i
        i_max = i
        for j in range(i+1,N):
          if (abs(A[i_max,i]) < abs(A[j,i])): i_max = j
        if (A[i_max,i] < 0.0):
            A[i_max,0:] = A[i_max,0:]*-1.0 # Switch signs
            A_inv[i_max,0:] = A_inv[i_max,0:]*-1.0
        # Move pivot row so column is on diagonal
        if (i_max != i): # Swap rows
            P[i_max, i_max] = 0.0
            P[i, i] = 0.0
            P[i_max, i] = 1.0 # Move row i to row i_max
            P[i, i_max] = 1.0 # Move row i_max to row i
            A = np.matmul(P,A)
            A_inv = np.matmul(P,A_inv) 
            P_gen = np.matmul(P, P_gen) # !compt. intensive!
        # Pivot
        fac = A[i,i]
        A_inv[i,0:] = A_inv[i,0:]/fac
        A[i,0:] = A[i,0:]/fac # pivot element becomes 1
        for j in range(i+1,N): # Other rows
            if  (A[j,i] != 0.0): # skip rows with column element = 0
                fac = A[j,i] # factor for multiplying pivot row
                A[j,0:] -= fac*A[i,0:]
                A_inv[j,0:] -= fac*A_inv[i,0:]
    return A 

rng2 = np.random.default_rng()
N = rng2.integers(6,14)
gen_num = rng2.integers(10,100000)
A, x, b = rand_SPD_Axb(N, gen_num) # Random SPD problem generator

A_red = row_reduc(A)

```

The random SPD coefficient matrix is:
<p align="center">
<img src="https://latex.codecogs.com/svg.image?A&space;=&space;\begin{bmatrix}11&space;&space;&&space;2&space;&space;&space;&&space;3.5&space;&&space;1.5&space;&&space;3.5&space;&&space;2&space;&space;&space;&&space;3&space;&space;&space;\\2&space;&space;&space;&&space;7&space;&space;&space;&&space;0.5&space;&&space;5&space;&space;&space;&&space;2.5&space;&&space;2&space;&space;&space;&&space;2.5&space;\\3.5&space;&&space;0.5&space;&&space;10&space;&space;&&space;2.5&space;&&space;3&space;&space;&space;&&space;4&space;&space;&space;&&space;1&space;&space;&space;\\1.5&space;&&space;5&space;&space;&space;&&space;2.5&space;&&space;8&space;&space;&space;&&space;0.5&space;&&space;2.5&space;&&space;3.5&space;\\3.5&space;&&space;2.5&space;&&space;3&space;&space;&space;&&space;0.5&space;&&space;12&space;&space;&&space;2&space;&space;&space;&&space;2&space;&space;&space;\\2&space;&space;&space;&&space;2&space;&space;&space;&&space;4&space;&space;&space;&&space;2.5&space;&&space;2&space;&space;&space;&&space;11&space;&space;&&space;1&space;&space;&space;\\3&space;&space;&space;&&space;2.5&space;&&space;1&space;&space;&space;&&space;3.5&space;&&space;2&space;&space;&space;&&space;1&space;&space;&space;&&space;12&space;\end{bmatrix}"/>
</p>
After passing it to the row_reduc function, we get:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?A&space;=&space;\begin{bmatrix}1&space;&&space;0.181818&space;&&space;0.318182&space;&&space;0.136364&space;&&space;0.318182&space;&&space;0.181818&space;&&space;0.272727&space;\\0&space;&&space;1&space;&space;&space;&space;&space;&space;&space;&space;&&space;-0.02055&space;&&space;0.712329&space;&&space;0.280822&space;&&space;0.246575&space;&&space;0.294521&space;\\0&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;1&space;&space;&space;&space;&space;&space;&space;&space;&&space;0.238628&space;&&space;0.216654&space;&&space;0.382421&space;&&space;0.009638&space;\\0&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;1&space;&space;&space;&space;&space;&space;&space;&space;&&space;-0.44976&space;&&space;0.063985&space;&&space;0.427869&space;\\0&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;1&space;&space;&space;&space;&space;&space;&space;&space;&&space;0.030697&space;&&space;0.134696&space;\\0&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;1&space;&space;&space;&space;&space;&space;&space;&space;&&space;-0.02305&space;\\0&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;0&space;&space;&space;&space;&space;&space;&space;&space;&&space;9.7165&space;&space;\end{bmatrix}"/>
</p>

The diagonal elements of the matrix in row-echelon form are its pivots. A matrix is symmetric positive-definite if those pivots are positive non-zero values, as we've confirmed above. This gives us options with how we want to solve the system.


## Solution Methods

There are a number of methods that can be applied to solving linear systems of equations with symmetric positive-definite coefficient matrices. We may classify these methods as iterative or direct. Iterative methods, like the Gauss-Seidel method we've already used, work by using successive solutions to get closer to the true solution. Apart from the Gauss-Seidel method, there is the similar Jacobi method as well as the more sophisticated Conjugate Gradient method. Direct methods solve the system Ax = b by finding the inverse of A and multiplying by the b matrix. Gaussian elimination, a variant of which is used in this post (the Gauss-Jordan method) is one method and Cholesky decomposition is another. 

The Jacobi method is like Gauss-Seidel but uses the the previous solution at each iteration. Dr. Zhiliang Xu from Notre Dame University presents a [brief comparison](https://www3.nd.edu/~zxu2/acms40390F12/Lec-7.3.pdf) of the two methods and discusses convergence rates.

I recommend to the reader L. Ridgeway Scott's discussion of the Conjugate Gradient method in their Numerical Analysis text [ref] (with Wikipedia for some orientation). Gross high-level simplification: imagine a space where the the true solution of a function Ax = b is separated from an initial conditon x_0. The conjugate gradient method finds a numerical solution close to the true solution by minimizing the distance between the numerical solution and the true solution (_Q<sub>A</sub> = x<sup>T</sup> Ax/2 - b<sup>T</sup>x_). To minimize this distance, its gradient _grad Q<sub>A</sub>_ and the numerical residual _Ax - b_ are used to find a direction and step size for for the numerical solution to be shifted. The [Wikipedia page](https://en.wikipedia.org/wiki/Conjugate_gradient_method) has a handy graphic showing this process geometrically and distinguishes the conjugate gradient from the gradient descent method (also seen in ML) by how the search directions differ: the conjugate gradient method changes each step's search direction so that sequential search directions are orthogonal (or conjugate) to each other. These direction changes avoid repeating solutions or search directions. 

In mathematical notation, we consider the problem _Cy=g_ and start with the residual as the search direction:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\begin{align*}Cy&space;&=&space;g&space;\\&space;s_0&space;=&space;r_0&space;&=&space;Cy_0-g\end{align*}" />
</p>

Starting with our initial condition _y<sub>0</sub>_, we can calculate the initial residual _r<sub>0</sub>_ which is set to the initial search direction _s<sub>0</sub>_. Then, the iteration begins:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\begin{align*}for\;\;k &= 0,1,2,...:\\ &\alpha_k =-\dfrac{(r_k,s_k)_I}{(s_k,z_k)_I},\;\; (r_k,s_k)_I = r_k^T s_k\\&y_{k+1} = y_k + \alpha_k s_k,\;\; z_k=Cs_k\\ &r_{k+1} = Cy_{k+1}-g=r_k+\alpha_kCs_k=r_k+\alpha_kz_k\\ &\beta_k = -\dfrac{(r_{k+1},z_k)_I}{(s_k,z_k)_I}\\ &s_{k+1} = r_{k+1} + \beta_k s_k \end{align*} " />
</p>

In the above algorithm, all lower-case latin letters are column vectors and greek letters are scalars. With each iteration, the algorithm moves the numerical solution _y<sub>k</sub>_ closer to the true solution by shifting it based on search directions _s<sub>k</sub>_, which itsef changes based on the residual and its previous iterations. The conjugate gradient method differs from other gradient descent methods by using search directions which are orthogonal to each other so that the solutions aren't repeated. 

Direct solution methods basically solve the system _Ax = b_ as _x = A<sup>-1</sup> b_, which requires calculating the inverse. Gauss-Jordan elimination is a direct method where equations are reduced into a form that allows substituion. For example, in the pressure equation for Darcy's law, the boundary conditions are given. Using Gauss-Jordan elimination, one would be able to reduce the system of equations to a form (see Upper Triangular form) where you may substitute the outlet pressure into an equation for the adjacent cell's pressure. You may then solve for that cell's pressure which will be used in the next cell's equation. 

The Gauss-Jordan elimination algorithm included in the code presented here includes full-pivoting, basically conducting row operations to transform the original matrix A into an identity matrix I. By conducting the same row operations on the original matrix A and an identity matrix I, one can find the inverse of matrix A which can then be used as a linear operator with the source/constant vector b (_x = A<sup>-1</sup>b_). If the terms in A are constant over the time of a transient simulation, then this operation only needs to be done once, at the beginning, which saves computational time.

Cholesky factorization is another direct method which divides the original coefficient matrix A into a lower triangular matrix and its conjugate transpose. It then conducts a two-step solution where an intermediate vector _y_ is solved for with the lower triangular matrix _L_ and the original source/constant vector _b_ (_Ly = b_), then the desired vector x is solved for using the conjugate transpose _L<sup>*</sup>_ and the intermediate vector (_L<sup>*</sup>x = y_). Wikipedia's and NumPy's/SciPy's pages on the Cholesky factorization method vis-a-vis solving linear equations are useful sources for how it works. 

## Solution Code

The Gauss-Seidel method and Gauss-Jordan methods are constructed in Python below. The Gauss-Seidel method, as it was coded in the previous post, could only solve 1D systems of equations, where only adjacent cells interacted with each other. However, the method can solve any system of equations where the coefficient matrix is symmetric-positive-definite. The Gauss-Seidel method code is updated to:

```
#New Gauss-Seidel Code

def gauss_seidel(A,N,b,x0,tol):
    # Gauss-Seidel: x0 is the initial solution
    # No inverse matrix to return, just iteration states

    ## Initialize
    x = x0
    dif = np.matmul(A,x)-b
    res = np.array([LA.norm(dif)],dtype=float)
    # tol = 1E-8
    iter = 0
    max_iter = 1E5
    # print('> Gauss-Seidel Initial Residual: {}'.format(res[0]))
    ## Iteration Loop
    res_len = len(res)
    while ((res[res_len-1] > tol) & (iter < max_iter)):
        iter += 1
        for i in range(0,N): # Loop through equations
            x_sum = b[i]
            for j in range(0,N): 
                if (i != j): x_sum -= A[i,j]*x[j]
            x[i] = x_sum/A[i,i]

        if ((iter % 10) == 0):
            dif = np.matmul(A,x)-b
            res_iter = np.array(LA.norm(dif))
            res = np.append(res,res_iter) # check this code
            res_len = len(res)
        # print('> Gauss-Seidel Residual at iter {}: {}'.format(iter,res_iter))

    if (res[res_len-1] <= tol): 
        print('>> Gauss-Seidel Converged at iteration',iter,\
            '! tol =',tol,', res =',res[res_len-1])
    else:
        print('> Gauss-Seidel Divergence at iteration',iter,\
            '! tol =',tol,', res =',res[res_len-1])
    print_sol('Gauss-Seidel_Sol',A,A,x,np.matmul(A,x),N)

    return x, res[res_len-1], iter
```

We can test it with a random symmetric-positive-definite coefficient matrix _A_ and random solution vector _x_. Multiplying _A_ and _x_, we get the constant vector _b_:

```
import numpy as np
.
.
.
def rand_SPD_Axb(N,gen_num):
    # Random A, x, and corresponding generation for Ax = b problem
    rng = np.random.default_rng(gen_num)
    # A: symmetric positive definite
    A_seed = rng.integers(low=-5, high=5, size=(N,N)).astype(np.float)
    A_seed_T = np.transpose(A_seed)
    A = (np.abs(A_seed) + np.abs(A_seed_T))*0.5
    A += np.diag(np.ones(N)*N)
    A_org = A
    print('Determinant of A:',LA.det(A))
    x = rng.integers(low=-N, high=N, size=N).astype(np.float)
    b = np.matmul(A,x)
    print_mat1('A_Org',A,N)
    print_mat3('Original',A,np.eye(N),np.eye(N),b,N)
    
    return A, x, b
.
.
.
rng2 = np.random.default_rng()
N = rng2.integers(6,14)
gen_num = rng2.integers(10,100000)
A, x, b = rand_SPD_Axb(N, gen_num)	
```	

In the current case, the matrix _A_ and vector _b_ are:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?A&space;=&space;\begin{bmatrix}10&space;&space;&&space;3&space;&space;&space;&&space;1.5&space;&&space;1.5&space;&&space;1.5&space;&&space;2.5&space;\\3&space;&space;&space;&&space;9&space;&space;&space;&&space;3&space;&space;&space;&&space;1.5&space;&&space;2.5&space;&&space;4.5&space;\\1.5&space;&&space;3&space;&space;&space;&&space;6&space;&space;&space;&&space;3&space;&space;&space;&&space;2&space;&space;&space;&&space;3&space;&space;&space;\\1.5&space;&&space;1.5&space;&&space;3&space;&space;&space;&&space;7&space;&space;&space;&&space;3.5&space;&&space;3.5&space;\\1.5&space;&&space;2.5&space;&&space;2&space;&space;&space;&&space;3.5&space;&&space;10&space;&space;&&space;0.5&space;\\2.5&space;&&space;4.5&space;&&space;3&space;&space;&space;&&space;3.5&space;&&space;0.5&space;&&space;10&space;\end{bmatrix},\;\;&space;b&space;=&space;\begin{bmatrix}-8.5&space;\\30&space;\\43.5&space;\\35&space;\\\11.5&space;\\\61&space;\end{bmatrix}"/>
</p>

We then run the Gauss-Seidel method for the system and check the L2-norm of the residual vector _r = ||Ax-b||_:

```
>>> x
array([-3.,  0.,  5.,  1.,  0.,  5.])
>>> x_gs
array([-3.00000000e+00,  1.35871561e-12,  5.00000000e+00,  1.00000000e+00,
       -5.07416331e-13,  5.00000000e+00])
>>> LA.norm(np.matmul(A,x_gs)-b) 
2.465053629363389e-12
>>> LA.norm(x_gs-x) 
2.4648988254742506e-12
```

The true and Gauss-Seidel numerical solutions look close, both by difference between the vectors and the L2-norm of hte difference between the solutions. Additionally, the residual is small enough where we can assume that the numerical solution has converged. We can check how the method is performs by looking at individual iterations:

![Gauss-Seidel Iterations](../../../assets/images/GS_Iter_Out.png "Gauss-Seidel Iterations")

In the line legend, the letters after _x_ are the iteration number. The solution, visually, is already very close after 10 iterations. In a real case, a closer solution can be reached by running the code longer, either by increasing the maximum number of iterations or lowering the tolerance (basically the residual required). 

### Gauss-Jordan

The Gauss-Jordan method is coded as:

```
# Gauss-Jordan Code

def gauss_jordan(A,N,b):
    # Custom Gauss-Jordan Elimination Algo

    # Notes from 9/11-13 2021 (Math folder)
    A_org = np.copy(A) # Original A (need for tests)
    b_org = np.copy(b)
    A_inv = np.eye(N) # A inverse initialized as identity matrix
    P_gen = np.eye(N) # General Permutation matrix
    P = np.eye(N) # Permutation matrix for this column pivot
    # Reduce to upper tridiagonal coefficient matrix
    for i in range(0, N-1):
        # Find row of pivot element in column i
        i_max = i
        for j in range(i+1,N):
          if (abs(A[i_max,i]) < abs(A[j,i])): i_max = j
        if (A[i_max,i] < 0.0):
            A[i_max,0:] = A[i_max,0:]*-1.0 # Switch signs
            A_inv[i_max,0:] = A_inv[i_max,0:]*-1.0
            b[i_max] = b[i_max]*-1.0
        # Move pivot row so column is on diagonal
        if (i_max != i): # Swap rows
            P[i_max, i_max] = 0.0
            P[i, i] = 0.0
            P[i_max, i] = 1.0 # Move row i to row i_max
            P[i, i_max] = 1.0 # Move row i_max to row i
            A = np.matmul(P,A)
            A_inv = np.matmul(P,A_inv) 
            b = np.matmul(P,b)
            P_gen = np.matmul(P, P_gen) # !compt. intensive!
        # Pivot
        fac = A[i,i]
        b[i] = b[i]/fac # A[i,i] is now pivot element
        A_inv[i,0:] = A_inv[i,0:]/fac
        A[i,0:] = A[i,0:]/fac # pivot element becomes 1
        for j in range(i+1,N): # Other rows
            if  (A[j,i] != 0.0): # skip rows with column element = 0
                fac = A[j,i] # factor for multiplying pivot row
                A[j,0:] -= fac*A[i,0:]
                A_inv[j,0:] -= fac*A_inv[i,0:]
                b[j] -= fac*b[i]
    if (A[N-1,N-1] != 0.0):
        fac = A[N-1,N-1]
        b[N-1] = b[N-1]/fac # last row
        A_inv[N-1,0:] = A_inv[N-1,0:]/fac
        A[N-1,0:] = A[N-1,0:]/fac
    elif (b[N-1] != 0.0): # No solution here
        raise Exception('Singular matrix: 0*x_N != 0.0')
    
    # Check: Matrix is in upper tridiagonal form
    if (LA.norm(A-np.triu(A)) <= 1.0E-14):
        print('>> A reduced to upper tridiagonal matrix! <<')
    else: 
        print('> Partial Reduction Failed!')
        return b, A
    
    # Back sub for trailing terms (full pivoting)
    for i in range(N-1,0,-1): # Start: bottom row cancels out other rows
        for j in range(i-1,-1,-1): # Other rows
            if (A[j,i] != 0.0):
                fac = A[j,i]
                A[j,0:] -= fac*A[i,0:]
                A_inv[j,0:] -= fac*A_inv[i,0:]
                b[j] -= fac*b[i]
    
    # Revert solutions back to original order
    P_trans = np.transpose(P)
    A_inv_perm = np.matmul(P_trans,A_inv)
   
    # Solution Test: norm(x - (A_inv x b)) < some precision number
    l2_res = LA.norm(b-np.matmul(A_inv,b_org))

    # Check: Full pivoting has turned A into identity matrix
    if (LA.norm(A-np.eye(N))<1.0E-14): 
        return b, A_inv
    else:
        print('> Full Pivoting Failed!')
        return b, A
```

So we run it and check results and residual:

```
>>> x
array([-3.,  0.,  5.,  1.,  0.,  5.])
>>> x_gj
array([-3.00000000e+00, -4.44089210e-16,  5.00000000e+00,  1.00000000e+00,
       -2.22044605e-16,  5.00000000e+00])
>>> LA.norm(np.matmul(A,x_gj)-b)  
0.0
>>> LA.norm(x_gj-x)                    
1.4895204919483639e-15
```

The Gauss-Jordan residual is zero and the L2-norm of the error vector is on the order of _1E-15_, which is orders of magnitude smaller than the Gauss-Seidel method. While the residuals are different, the solutions will appear to overlap.  For the purposes of this test, they are close enough. However, for a transient study, where a series of solutions will be stacked on top of each other, this residuals become important as the errors will accumulate. 

So far, we've coded the schemes explicitly, which forces us to understand the schemes but also allows us to report more data from the solver, for example, how the solver behaves at intermediate solutions. This can be invaluable if we're not sure whether our understanding of the problem is accurate. However, once we're sure that we understand the mathematical problem and that the coefficient matrix is correct, we can use solvers already in the Python/NumPy/SciPy libraries to get the solution. The advantages, beyond not having to debug our own code, are the optimizations that the package's developers have implemented. Though not having to debug is a huge advantage when solving small toy problems, the optimizations can be crucial when dealing with multi-billion calculation problems.

We call the conjugate gradient method from Sparse Linear Algebra toolbox from SciPy and the Cholesky factorization method from SciPy's Linear Algebra toolbox as:

```
import scipy.sparse.linalg as SLA
from scipy.linalg import cho_factor, cho_solve
.
.
.

[x_cg, info_cg] = cg(A, b)    
c_cho, low = cho_factor(A)
x_cho = cho_solve((c_cho, low), b)
```

The results for the conjugate gradient method are:

```
>>> x_cg
array([-3.00000000e+00, -4.44089210e-16,  5.00000000e+00,  1.00000000e+00,
       -2.22044605e-16,  5.00000000e+00])
>>> LA.norm(np.matmul(A,x_cg)-b)
0.0
>>> LA.norm(x_cg-x)              
1.4895204919483639e-15
```

and, for the Cholesky factorization:

```
>>> x_cho
array([-3.00000000e+00, -4.44089210e-16,  5.00000000e+00,  1.00000000e+00,
       -2.22044605e-16,  5.00000000e+00])
>>> LA.norm(np.matmul(A,x_cho)-b)
0.0
>>> LA.norm(x_cho-x)              
1.4895204919483639e-15
```

The results show about the same accuracy as the Gauss-Jordan method we coded, but without the having to code.


## Apply to Darcy\'s Law

Now that we know our four methods work on the general class of problems (symmetric-positive-definite coefficient matrix), we can now apply them to our specific problem and see which one works best. We set up the terms for the pressures in the [previous post](https://tim-munuhe.github.io/2021/06/07/solve-1d-darcy-continuity). This post will leave it to the reader to put the right terms for the pressure equations from the previous post into the coefficient matrix at the beginning of this post. Just know, the problem has not changed, just the form. The python code may also be deciphered for clues as to the terms.

We call methods on the coefficient matrices and constant vectors:

```
import multi_solv_v3 as ms
import darcy_v1 as dc
import scipy.sparse.linalg as SLA
from scipy.linalg import cho_factor, cho_solve
import csv
import numpy as np
from numpy import linalg as LA
import copy
import time

# Read Case
with open('casefile_mesh.csv',newline='') as casefile:  # open and read case params
    casereader = csv.DictReader(casefile) 
    i = 0 
    caselist = {} 
    for row in casereader:    
        caselist[i] = row 
        i += 1         

for i in range(0,len(caselist)):
      # Create case object and check values of variables
      case_current = dc.case_param(caselist[i])
      print(case_current.x0)
      print(case_current.dx)

      # Initialize and check mesh object creation ##
      case_mesh = dc.mesh(case_current) # create case mesh from case parameters
      Nx = case_mesh.Nx
      print('Node Locations w/ inlet:', case_mesh.x[0:5]) # check inlet location and spacing
      print('Nx:', case_mesh.Nx) # check number of elements
      print('Outlet Location:', case_mesh.x[case_mesh.Nx-1])
      print('Face Locations:', case_mesh.xc[0:5]) # 

      # Create fluid and porous medium objects for this specific case ##
      fl1 = dc.fluid(case_mesh,case_current.fl) # fluid object, determined by mesh and case's fluid properties
      pm1 = dc.por_med(case_mesh,case_current.pm) # porous medium object, determined by mesh and case's porous medium properties

      # Linear P
      print('Original P:',fl1.p[0:4],fl1.p[case_mesh.Nx-5:case_mesh.Nx-1])
      fl1.p_lin(case_mesh)
      print('Linear P:',fl1.p[0:4],fl1.p[case_mesh.Nx-5:case_mesh.Nx-1])

      # Darcy Velocity
      print('Original u:',fl1.u[0:4]) # velocity from initialization
      fl1.u = np.zeros(case_mesh.Nx) # zero out velocity 
      fl1.darcyv(case_mesh,pm1) # use darcyv method
      print('Darcy u:',fl1.u[0:4]) # print to confirm that darcyv did what it was supposed to (got same solution as initialization)

      ## Pressure calculation using Different Solvers

      # Gauss-Seidel Original
      fl_gso = dc.fluid(case_mesh,case_current.fl)
      print('Original P:',fl_gso.p[0:4],
            fl_gso.p[Nx-5:Nx-1])
      A,b = fl_gso.coeff_Ab(case_mesh,pm1)
      time_gso = time.perf_counter()
      [iter_gso, res_gso] = fl_gso.gauss_seidel(case_mesh,pm1,1.0E-9)
      time_gso -= time.perf_counter()
      time_gso = abs(time_gso)
      fl_gso.darcyv(case_mesh,pm1)

      # Gauss-Seidel
      fl_gs = dc.fluid(case_mesh,case_current.fl)
      # print('Original P:',fl_gs.p[0:4],
      #       fl_gs.p[Nx-5:Nx-1])
      A,b = fl_gs.coeff_Ab(case_mesh,pm1)
      time_gs = time.perf_counter()
      fl_gs.p, res_gs, iter_gs = ms.gauss_seidel(A,Nx,b,np.zeros(Nx),1.0E-9)
      time_gs -= time.perf_counter()
      time_gs = abs(time_gs)
      fl_gs.darcyv(case_mesh,pm1)

      # Gauss-Jordan
      fl_gj = dc.fluid(case_mesh,case_current.fl)
      time_gj = time.perf_counter()
      fl_gj.p, res_gj = ms.gauss_jordan(A,Nx,b)
      time_gj -= time.perf_counter()
      time_gj = abs(time_gj)
      fl_gj.darcyv(case_mesh,pm1)

      # Conjugate Gradient
      fl_cg = dc.fluid(case_mesh,case_current.fl)
      time_cg = time.perf_counter()
      fl_cg.p, info_cg = SLA.cg(A, b, tol=1.0E-9)
      time_cg -= time.perf_counter()
      time_cg = abs(time_cg)
      fl_cg.darcyv(case_mesh,pm1) 

      # Conjugate Gradient Squared
      fl_cgs = dc.fluid(case_mesh,case_current.fl)
      time_cgs = time.perf_counter()
      fl_cgs.p, info_cg = SLA.cgs(A, b, tol=1.0E-9)
      time_cgs -= time.perf_counter()
      time_cgs = abs(time_cgs)
      fl_cgs.darcyv(case_mesh,pm1) 

      # Cholesky Factorization
      fl_cho = dc.fluid(case_mesh,case_current.fl)
      time_cho = time.perf_counter()
      c_cho, low = cho_factor(A)
      fl_cho.p = cho_solve((c_cho, low), b)
      time_cho -= time.perf_counter()
      time_cho = abs(time_cho)
      fl_cho.darcyv(case_mesh,pm1)

```


The L2-norms of the error vectors show how poor the Gauss-Seidel method compares against the other methods:

```
>>> LA.norm(fl_gs.p-fl1.p)
0.01010884570852091
>>> LA.norm(fl_gj.p-fl1.p)  
5.312732304408523e-12
>>> LA.norm(fl_cg.p-fl1.p)  
5.312732304408523e-12
>>> LA.norm(fl_cho.p-fl1.p)  
5.312732304408523e-12
```

Apart from the Gauss-Seidel method, the methods have reasonably close residuals. What does this look like with regards to the pressure fields?

![P vs. x for different methods.](../../../assets/images/p_N100.png "P vs. x for different methods.")

They virtually overlap, despite the large difference in the error norms. So, apart from the error norm, how else would we distinguish the methods? Since our case is relatively simple, we can use two basic criteria: timing and memory usage. Both are important for very large problems, for example if we broke the domain into millions of node elements. In comparison, the tests we did earlier with the random symmetric-positive-definite coefficient matrix only used 6 elements. For this post, only time will be considered:

```
Timing (s): G-Seidel | G-Jordan | Conj. Grad | Cholesky
----------------------------------------------------------------
79.4208762 | 0.0295643 | 0.0014511 | 0.0019525
```

With _100_ elements, it takes over a minute for the Gauss-Seidel method while the other methods take less than one second. 

So why would we use the Gauss-Seidel method? Using the method as coded in the previous post, with the equations coded directly (center, upstream, and downstream terms  rather than put into _Ax = b_ form, the only values kept in memory were the pressure vector, permeability, and viscosity, a total of _3N_ terms. The other methods have a minimum of _N<sup>2</sup>_ values (assuming sparse matrices aren't used). In the case that the domain is divided into millions of elements, it may only be possible to use the Gauss-Seidel method to solve the problem on a personal computer. 