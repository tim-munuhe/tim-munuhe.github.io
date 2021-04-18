---
layout: post
title: Basic Python with Darcy's Law
description: Darcy's Law for Porous Fluid Flow using Python
summary: Darcy's Law is a simple enough equation to use with Python and is a useful pier for exploring scientific computing.
tags: [css]
---

DRAFT/OUTLINE

-	Darcy’s Law: used to calculate flow velocity in porous medium under a pressure gradient
-	Useful test problem to solve because it is PDE, class of equations with wide variety of uses
	o	Other famous PDEs (Wave equation, Heat equation, Black-Scholes equation, etc.
-	Darcy’s Law allows physical intuition and can be simplified or complicated with relatively intuitive variables
-	Can add to be part of greater model (see later)
-	Problem Statement: Horizontal Channel filled with sand (no gravity)

-> Darcy's Law here; u = -K/mu gradP (Eq. 1)

```css
.♫ { background: #222; color: #FFF; } 
.ಠ_ಠ { background: #ccc; color: #fff; }
```

-	If we know the pressure field and sand properties, we can solve this pretty simply
	o	Assume the sand is homogeneous and that water is flowing through the pores: Can find the K and mu with various relations and tables
	o	If the sand is homogeneous and the water does not change properties, we can also assume the pressure gradient is constant (come back to that later)
	o	We just need the pressure conditions at the inlet and outlet

-> Pressure distro here (Eq. 2)
-> Pressure gradient here (Eq. 3)
-> superficial velocity u = -K/mu*(P_L-P_0)/L (Eq. 4)


-	Equation 4 is the superficial velocity. In the simplified problem that we’re considering, we can multiply this by the area of the inlet/outlet to get the flow rate (what folks actually measure). 

```css
a { text-decoration: overline red wavy; }
```

Solving for velocity over a mesh with known pressure gradient

-	It’s unlikely you need python for this. More likely, you’ll need to code more complex problems that you can’t do (efficiently) with pen and paper: heterogeneous porous medium
	Solve for the velocity at every point within the channel: you can see where flow is fastest/slowest depending on the porous medium parameters
-	One step: you still know the pressure distribution throughout the channel but want to see how velocity changes if you change the pressure distribution
	Let’s create a mesh to represent the relevant properties at specified positions throughout the porous medium
	o	Use Finite Difference method, nodes and faces
	o	Staggered mesh: nodes have pressure, faces have velocity
	Check to see that the solution from the analytical part is the same: 
	o	use equation 2 for pressure at nodes then calculate velocity at faces
	o	use equation 1 for velocity at faces  discretize that pressure gradient term (math)
	o	check pressure gradient with equation 3
	Test it with a new unphysical pressure distribution and see pressure gradient and velocity

