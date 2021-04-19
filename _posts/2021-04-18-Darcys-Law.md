Basic Python with Darcy's Law
-----------------------------

Darcy's law is an equation used to calculate fluid flow through a porous medium under a pressure gradient. It's also a useful test problem to explore scientific computing with Python from a basic to intermediate level.

Darcy's Law can be written as: \\[\\vec{u} = -\\frac{K}{\\mu} \\nabla P \\]

where \\(\\vec{u}\\) is the superficial flow velocity, \\(K\\) is the hydraulic permeability of the porous medium, \\(\\mu\\) is the viscosity of the liquid, and \\(P\\) is the pressure distribution. All three variables can vary over space. If the problem considered is 2D or 3D, then Darcy's Law becomes a partial differential equation (PDE).

PDEs represent a myriad of phenomena mathematically, including heat transfer, electro-magnetism and the price of European options. You can find better discussions elsewhere but some prototypical equations to learn about are Laplace's equation, Poisson's equation, the Heat equation, and the Wave equation.

### Problem Setup

Let's say we have a pipe filled with sand. A fluid can flow in the spaces between the individual grains of sand, termed the pores. Depending on how big or small or well-packed the sand grains are, it's easier or harder for the fluid to flow through the pipe. This is represented by \\(K\\). The fluid's viscosity also affects how easily it can flow through the pores (maple syrup? or water? or air?).

Basically, the porous medium and fluid flow properties are constant. Let's also assume that the pipe is long enough relative to its diameter that we can assume that pressure only varies significantly along its axis. Then, we can treat this as a 1D problem: \\[u = -\\frac{K}{\\mu} \\dfrac{dP}{dx} \\]

So now we have an ordinary differential equation, or ODE. To complete the description of the problem we need 2 boundary conditions. Let's give two boundary conditions: \\[P(x=0) = P\_0\\] \\[P(x=L) = P\_L\\]

Let's give some properties so we can move on:

||
|\\(K\\)|\\(1\\times 10\^{-9} \\dfrac{m\^2}{s}\\)|
|\\(\\mu\\)|\\(0.001 Pa\\cdot s\\)|
|\\(P\_0\\)|\\(0 Pa\\)|
|\\(P\_L\\)|\\(-100 Pa\\)|
|\\(L\\)|\\(1 m\\)|

### Enough Physics. Let's Code!

With our current assumptions, the superficial velocity at every point within the pipe is: \\[u = -\\dfrac{K}{\\mu} \\dfrac{P\_L - P\_0}{L} \\] Then we can solve it quite easily with:

     K = 1.0E-9 # permeability  mu = 0.001 # viscosity  P_0 = 0.0 # inlet  P_L = -100.0 #outlet  L = 1.0 # pipe length  u = -K/mu*(P_L-P_0)/L  print(u)  >>> 9.999999999999999e-05   

Simple, but now I can change the variables and get the velocity immediately. I can even add a bit of extra code to output results and create a sort of solution space examining the effects of different variables. However, the assumptions made to get here are pretty restrictive. What if we want to check the pressure along the pipe? What if the sand is not homogeneous? What if the viscosity of the fluid changes because of temperature? We need a more robust solution.

### Some basic object-oriented programming

Admittedly, I'm still learning object-oriented programming so all I can do is write how I understand my code within the paradigm. That being said, the code will work, so take solace in that.

I want a more robust code that can take user input and tell the rest of the code how to run. So, I'm going to create a case object. I create a case class:

    class case_param():     def __init__(self):         self.dim = 1 # dimensions         self.x0 = 0.0 # inlet position         self.xL = 1.0 # outlet         fluid_name = 'Water'         mu = 0.001         u0 = 0.0         p0 = 0.0 # inlet pressure         pL = -100.0 # outlet         self.fl = {'Name': fluid_name, 'mu': mu, 'u0': u0, 'p0': p0, 'pL': pL}         pm_name = 'Sand'         K = 1.0E-9         eps = 0.15         self.pm = {'Name': pm_name, 'K':K, 'eps':eps}  

This class, in short, defines case objects through common variables: number of dimensions, inlet and outlet position, and the fluid and porous medium used and their properties. The fluid and porous medium are both represented thorugh dictionaries which other objects or methods can refer to. Let's use it:

    base = case_param() base.u0 = -base.pm['K']/base.fl['mu']*(base.fl['pL']-base.fl['p0'])/(base.xL-base.x0) print(base.u0) >>> 9.999999999999999e-05  

So, we've got the same result as the previous, simpler code. We've also created a case object that the other to-be-created code can use. For now, let me compress the code by initializing the superficial velocity in the instantiation (\_\_init\_\_ method):

    class case_param():     def __init__(self):         self.dim = 1 # dimensions         self.x0 = 0.0 # inlet position         self.xL = 1.0 # outlet         fluid_name = 'Water'         mu = 0.001         u0 = 0.0         p0 = 0.0 # inlet pressure         pL = -100.0 # outlet         self.fl = {'Name': fluid_name, 'mu': mu, 'u0': u0, 'p0': p0, 'pL': pL}         pm_name = 'Sand'         K = 1.0E-9         eps = 0.15         self.pm = {'Name': pm_name, 'K':K, 'eps':eps}         self.fl['u0'] = -K/mu*(pL-p0)/(self.xL-self.x0) base = case_param() print(base.fl['u0']) >>> 9.999999999999999e-05  

Later, I'll be able to pass this to a mesh object and Darcy's law method to create a solution that I can plot and output to CSV.

April 18, 2021 ·
