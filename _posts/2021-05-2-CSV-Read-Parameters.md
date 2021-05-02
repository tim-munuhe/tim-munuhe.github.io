---
layout: post
title: Basic Python with Darcy's Law: User-Defined Parameters
description: Read parameters for Darcy's Law from a CSV file.
summary: Make our future work more streamlined by reading multiple case's parameters from a CSV file.
mathjax: true
---

In the previous blog post, I introduced some object-oriented Python by creating a parameter class that initialized a solution to Darcy\'s law for some default fluid and porous medium. So, just by creating the case, we were done. Now, we want to feed our own properties into Darcy\'s Law. 

Let\'s start with the (easier) parameter problem: we can feed in our own fluid and porous medium properties by adding to the instantiation's input. Currently, we just have (self) for the initialization. Let\'s add the length of the porous medium:

<pre><code class="language-python">class case_param(): <br>    def __init__(self,L): <br>        self.dim = 1 # dimensions <br>        self.x0 = 0.0 # inlet position <br>        self.xL = self.x0 + L # outlet <br>        fluid_name = 'Water' <br>        mu = 0.001 <br>        u0 = 0.0 <br>        p0 = 0.0 # inlet pressure <br>        pL = -100.0 # outlet <br>        self.fl = {'Name': fluid_name, 'mu': mu, 'u0': u0, 'p0': p0, 'pL': pL} <br>        pm_name = 'Sand' <br>        K = 1.0E-9 <br>        eps = 0.15 <br>        self.pm = {'Name': pm_name, 'K':K, 'eps':eps} >br>         self.fl['u0'] = -K/mu*(pL-p0)/(self.xL-self.x0) </code> </pre>

When we create the case_param object, we need to give the length in the parentheses, like so:

<pre><code class="language-python">base = case_param(1.0)

Then, when we call the case\'s outlet location, we get our length:

<pre><code class="language-python">print(base.xL) <br> >>> 9.999999999999999e-05  </code> </pre>

which gives us the same answer as the first blog post:

<pre><code class="language-python">print(base.fl['u0']) <br> >>> 9.999999999999999e-05  </code> </pre>

We could keep going this way but we have at least 8 parameters for our case that we\'d want to vary: the fluid name, viscosity, inlet and outlet pressure, porous medium name, permeability, porosity, and the length of the domain. Instead, we can use text files or CSV files with specified formats to feed in the case parameters, allowing a more streamlined multi-case process. Let\'s use CSV: we can create it in Excel and it's use in Python for Data Science means there will be resources for troublshooting later (;-)). 

First, import the CSV package:

<pre><code class="language-python">print(base.fl['u0']) <br> >>> 9.999999999999999e-05  </code> </pre>

Next, let's create our CSV case file using Excel:

| case\_name | fluid | p0      | pL        | mu    | porous\_medium | length | K        | eps   |
| ---------- | ----- | ------- | --------- | ----- | -------------- | ------ | -------- | ----- |
| base       | water | 0.000   | \-100.000 | 0.001 | sand           | 1.000  | 1.00E-09 | 0.150 |
| long       | water | 0.000   | \-100.000 | 0.001 | sand           | 2.000  | 1.00E-09 | 0.150 |
| press      | water | 100.000 | \-100.000 | 0.001 | sand           | 1.000  | 1.00E-09 | 0.150 |
| powder     | water | 0.000   | \-100.000 | 0.001 | powder         | 1.000  | 1.00E-11 | 0.300 |
| oil        | oil   | 0.000   | \-100.000 | 0.060 | sane           | 1.000  | 1.00E-09 | 0.150 |

We can use the csv.reader function and skip the first line to create individual case parameter lists, or, we can use the csv.DictReader function to construct individual case dictionaries:

<pre><code class="language-python">with open('casefile.csv',newline='') as casefile: <br>    casereader = csv.DictReader(casefile) <br>    i = 0 <br>    caselist = {} <br>    for row in casereader:    <br>        caselist[i] = row <br>        print(row['case_name'], row['fluid'], row['mu']) # check that code works as expected <br>        i += 1 <br> >>> base water 0.001 <br> >>> long water 0.001 <br> >>> press water 0.001 <br> >>> powder water 0.001 <br> >>> oil oil 0.060 </code> </pre>

DictReader uses the first row of the CSV file as the keys and the subsequent row values are the dictionary entries. The only problem is that all entries are read as strings, which must convert the number variables to floats in the <pre><code class="language-python">case_param</code> </pre> instantiation:

<pre><code class="language-python">class case_param(): <br>    def __init__(self,param): <br>        self.name = param['case_name'] # now the name is given inside the case, not as the case's actually name <br>        self.dim = 1 # dimensions <br>        self.x0 = 0.0 # inlet position <br>        self.xL = self.x0 + float(param['length']) # outlet <br>        fluid_name = param['fluid'] <br>        mu = float(param['mu']) <br>        u0 = 0.0 <br>        p0 = float(param['p0']) # inlet pressure <br>        pL = float(param['pL']) # outlet <br>        self.fl = {'Name': fluid_name, 'mu': mu, 'u0': u0, 'p0': p0, 'pL': pL} <br>        pm_name = param['porous_medium']  <br>        K = float(param['K']) <br>        eps = float(param['eps']) <br>        self.pm = {'Name': pm_name, 'K':K, 'eps':eps} <br>         self.fl['u0'] = -K/mu*(pL-p0)/(self.xL-self.x0) </code> </pre>

We can initialize the original base case and the oil case and compare the velocities:

<pre><code class="language-python"> base = case_param(caselist[0]) <br> oil = case_param(caselist[4]) <br> print(base.fl['u0'] <br> >>> 9.999999999999999e-05 <br> print(oil.fl['u0'] <br> >>> 1.6666666666666667e-06 </code> </pre>

We can see that the viscous oil slows down the flow, as expected. We haven't changed the Darcy\'s Law calculation so we know, at least, that the code is reading the CSV file correctly and initializing the case properly. 

The next step is to see what is going on between the inlet and outlet, specifically with the pressure.