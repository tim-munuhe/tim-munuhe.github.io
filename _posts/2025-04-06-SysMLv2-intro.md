---
layout: post
title: SysML v2 w/ Jupyter Hub - Intro
description: Let's take a look at SysML v2 model development and visualization using Jupyter Hub
summary: SysML v2 introduction
mathjax: true
---

SysML v2 is...

Math placeholder (Latex):
<p align="center">
<img src="https://latex.codecogs.com/svg.image?\vec{u}&space;=&space;-\frac{K}{\mu}\nabla&space;P"/>
</p>

where ![vel](https://latex.codecogs.com/svg.image?\vec{u}) is the superficial flow velocity, ![perm](https://latex.codecogs.com/svg.image?K) is the hydraulic permeability of the porous medium, ![visc](https://latex.codecogs.com/svg.image?\mu) is the viscosity of the liquid, and ![pres](https://latex.codecogs.com/svg.image?P) is the pressure distribution. All three variables can vary over space. If the problem considered is 2D or 3D, then Darcy\'s Law becomes a partial differential equation (PDE).



<h3>Heading3</h3>

Text

Let\'s give some properties so we can move on:

Property | Value
:-: | :-: 
![perm](https://latex.codecogs.com/svg.image?K) | ![Kval](https://latex.codecogs.com/svg.image?10^{-9}\;\frac{m^2}{s})
![visc](https://latex.codecogs.com/svg.image?\mu) | ![Kval](https://latex.codecogs.com/svg.image?0.001\;Pa\cdot&space;s)
![perm](https://latex.codecogs.com/svg.image?P_0) | ![Kval](https://latex.codecogs.com/svg.image?0\;Pa)
![perm](https://latex.codecogs.com/svg.image?P_L) | ![Kval](https://latex.codecogs.com/svg.image?-100\;Pa)
![perm](https://latex.codecogs.com/svg.image?L) | ![Kval](https://latex.codecogs.com/svg.image?1\;m)



<h3>Heading3</h3>

With our current assumptions, the superficial velocity at every point within the pipe is: 

<p align="center">
<img src="https://latex.codecogs.com/svg.image?u=-\dfrac{K}{\mu}\dfrac{P_L-P_0}{L}"/>
</p> 

Code placeholder:

<pre><code class="language-python"> K = 1.0E-9 # permeability <br> mu = 0.001 # viscosity <br> P_0 = 0.0 # inlet <br> P_L = -100.0 #outlet <br> L = 1.0 # pipe length <br> u = -K/mu*(P_L-P_0)/L <br> print(u) <br> >>> 9.999999999999999e-05  </code> </pre>

<h3>Heading3</h3>
Class code pleaceholder:
<pre><code class="language-python">class case_param(): <br>    def __init__(self): <br>        self.dim = 1 # dimensions <br>        self.x0 = 0.0 # inlet position <br>        self.xL = 1.0 # outlet <br>        fluid_name = 'Water' <br>        mu = 0.001 <br>        u0 = 0.0 <br>        p0 = 0.0 # inlet pressure <br>        pL = -100.0 # outlet <br>        self.fl = {'Name': fluid_name, 'mu': mu, 'u0': u0, 'p0': p0, 'pL': pL} <br>        pm_name = 'Sand' <br>        K = 1.0E-9 <br>        eps = 0.15 <br>        self.pm = {'Name': pm_name, 'K':K, 'eps':eps} </code> </pre>

Class Usage Code Placeholder
<pre><code class="language-python">base = case_param() <br>base.u0 = -base.pm['K']/base.fl['mu']*(base.fl['pL']-base.fl['p0'])/(base.xL-base.x0) <br>print(base.u0) <br>>>> 9.999999999999999e-05 </code> </pre>

More code: (\__init__ method):
<pre><code class="language-python">class case_param(): <br>    def __init__(self): <br>        self.dim = 1 # dimensions <br>        self.x0 = 0.0 # inlet position <br>        self.xL = 1.0 # outlet <br>        fluid_name = 'Water' <br>        mu = 0.001 <br>        u0 = 0.0 <br>        p0 = 0.0 # inlet pressure <br>        pL = -100.0 # outlet <br>        self.fl = {'Name': fluid_name, 'mu': mu, 'u0': u0, 'p0': p0, 'pL': pL} <br>        pm_name = 'Sand' <br>        K = 1.0E-9 <br>        eps = 0.15 <br>        self.pm = {'Name': pm_name, 'K':K, 'eps':eps} <br>        self.fl['u0'] = -K/mu*(pL-p0)/(self.xL-self.x0) <br><br>base = case_param() <br>print(base.fl['u0']) <br>>>> 9.999999999999999e-05 </code> </pre>
Later, I\'ll be able to pass this to a mesh object and Darcy's law method to create a solution that I can plot and output to CSV.
