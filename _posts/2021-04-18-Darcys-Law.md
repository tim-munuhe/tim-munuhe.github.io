<meta charset="UTF-8">

<meta name="viewport" content="width=device-width, initial-scale=1.0">

<meta name="generator" content="Jekyll v4.0.0">

<meta property="og:title" content="Basic Python w/ Darcy's Law">

<meta property="og:locale" content="en_US">

<meta name="description" content="This is a collection of short CSS snippets I thought might be useful for beginners">

<meta property="og:description" content="This is a collection of short CSS snippets I thought might be useful for beginners">

<link rel="canonical" href="http://localhost:4000/2014/05/12/css-hacks-you-may-not-know">

<meta property="og:url" content="http://localhost:4000/2014/05/12/css-hacks-you-may-not-know">

<meta property="og:site_name" content="Sidey">

<meta property="og:type" content="article">

<meta property="article:published_time" content="2014-05-12T00:00:00+03:00">

<meta name="twitter:card" content="summary">

<meta property="twitter:title" content="Basic Python w/ Darcy's Law">

<meta name="twitter:site" content="@">

<title> Basic Python w/ Darcy's Law </title>

<link rel="shortcut icon" href="/favicon.png">

<link rel="alternate" type="application/atom+xml" title="Sidey" href="/atom.xml">

<link rel="alternate" type="application/json" title="Sidey" href="http://localhost:4000/feed.json">

<link rel="sitemap" type="application/xml" title="sitemap" href="/sitemap.xml">

<style> *,:after,:before{box-sizing:border-box;background-color:inherit;color:inherit;margin:0;padding:0}body{font-family:-apple-system,BlinkMacSystemFont,'avenir next',avenir,helvetica,'helvetica neue',ubuntu,roboto,noto,'segoe ui',arial,sans-serif;-webkit-font-smoothing:antialiased;text-rendering:optimizeLegibility;line-height:1.5;font-size:1rem;color:#16171a}nav ul{border-right:1px solid #edf2f7}a{color:#000;text-decoration-skip-ink:auto;text-decoration:underline}pre{margin:.5rem 0;padding:.5rem}.post p{margin:.5rem 0}.post h1,.post h2,.post h3,.post h4{margin:1rem 0}.post h2:first-child,.project h2:first-child,.photo h2:first-child{margin-top:0}.meta{margin:2rem 0}code,pre{background:#ecedee}code{padding:.1rem}pre code{border:none}pre{padding:1rem;overflow-x:auto}img{max-width:100%}hr{background:#000;height:1px;border:0}header{flex-basis:10rem;flex-grow:1;position:relative}header a{text-decoration:none}header li{margin-bottom:.2rem;text-align:right;margin-right:2rem}header a.active{font-weight:bold}header,section{padding:1rem}blockquote{font-style:italic;border-left:5px solid #ececec;padding-left:1rem}h1,h2,h3,h4,h5{line-height:1;margin:1rem 0;font-weight:600}section h1:first-child{margin-top:0}strong,b{font-weight:bold}.photos ul{list-style:none}.photos li{margin-bottom:1.5rem}.photo picture,.project picture{margin-bottom:0.5rem}.posts ul,header ul{list-style:none}.posts li{align-items:center;display:flex;justify-content:space-between;margin-bottom:.5rem}.posts li a,.posts li div,.projects li a{white-space:nowrap;overflow:hidden;text-overflow:ellipsis;text-decoration:none}.posts li time,.projects li time{padding-left:1rem;white-space:nowrap;font-variant-numeric:tabular-nums}.post ul,.project ul,.post ol{list-style-position:inside}main{display:flex;flex-wrap:wrap;max-width:60rem;margin:2rem auto;padding:1rem}@media screen and (max-width: 45rem){header li{display:inline;margin-right:1rem}.logo{padding-bottom:1rem}header ul{border-bottom:1px solid #edf2f7;padding-bottom:2rem}nav ul{border-right:0px}.photos ul{margin-top:0.5rem}}section{flex-basis:0;flex-grow:999;min-width:70%;display:flex;flex-direction:column} table, th, td { border: 1px solid black; border-collapse: collapse; padding: 15px; text-align: center; } table.center { margin-left: auto; margin-right: auto; } </style>

<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>

<script id="MathJax-script" async="" src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

<main role="main"><header role="banner"><!--<h1 class="logo">Sidey</h1>--><nav role="navigation"><ul><li><a href="/">Writing</a></li><li><a href="/about">About</a></li><li><a href="/search">Search</a></li><li><a href="/atom.xml">Rss</a></li></ul></nav></header><section class="post"><h2>Basic Python with Darcy's Law </h2><p> Darcy's law is an equation used to calculate fluid flow through a porous medium under a pressure gradient. It's also a useful test problem to explore scientific computing with Python from a basic to intermediate level. </p><p> Darcy's Law can be written as: \[\vec{u} = -\frac{K}{\mu} \nabla P \] </p><p>where \(\vec{u}\) is the superficial flow velocity, \(K\) is the hydraulic permeability of the porous medium, \(\mu\) is the viscosity of the liquid, and \(P\) is the pressure distribution. All three variables can vary over space. If the problem considered is 2D or 3D, then Darcy's Law becomes a partial differential equation (PDE).</p><p> PDEs represent a myriad of phenomena mathematically, including heat transfer, electro-magnetism and the price of European options. You can find better discussions elsewhere but some prototypical equations to learn about are Laplace's equation, Poisson's equation, the Heat equation, and the Wave equation. </p><h3>Problem Setup</h3><p>Let's say we have a pipe filled with sand. A fluid can flow in the spaces between the individual grains of sand, termed the pores. Depending on how big or small or well-packed the sand grains are, it's easier or harder for the fluid to flow through the pipe. This is represented by \(K\). The fluid's viscosity also affects how easily it can flow through the pores (maple syrup? or water? or air?).</p><p>Basically, the porous medium and fluid flow properties are constant. Let's also assume that the pipe is long enough relative to its diameter that we can assume that pressure only varies significantly along its axis. Then, we can treat this as a 1D problem: \[u = -\frac{K}{\mu} \dfrac{dP}{dx} \]</p><p> So now we have an ordinary differential equation, or ODE. To complete the description of the problem we need 2 boundary conditions. Let's give two boundary conditions: \[P(x=0) = P_0\] \[P(x=L) = P_L\] </p><p>Let's give some properties so we can move on:</p><table class="center"><tbody><tr><td>\(K\)</td><td>\(1\times 10^{-9} \dfrac{m^2}{s}\)</td></tr><tr><td>\(\mu\)</td><td>\(0.001 Pa\cdot s\)</td></tr><tr><td>\(P_0\)</td><td>\(0 Pa\)</td></tr><tr><td>\(P_L\)</td><td>\(-100 Pa\)</td></tr><tr><td>\(L\)</td><td>\(1 m\)</td></tr></tbody></table><h3>Enough Physics. Let's Code!</h3><p>With our current assumptions, the superficial velocity at every point within the pipe is: \[u = -\dfrac{K}{\mu} \dfrac{P_L - P_0}{L} \] Then we can solve it quite easily with:</p><precode language="python" precodenum="0"></precode><p>Simple, but now I can change the variables and get the velocity immediately. I can even add a bit of extra code to output results and create a sort of solution space examining the effects of different variables. However, the assumptions made to get here are pretty restrictive. What if we want to check the pressure along the pipe? What if the sand is not homogeneous? What if the viscosity of the fluid changes because of temperature? We need a more robust solution. </p><h3>Some basic object-oriented programming</h3><p> Admittedly, I'm still learning object-oriented programming so all I can do is write how I understand my code within the paradigm. That being said, the code will work, so take solace in that. </p><p> I want a more robust code that can take user input and tell the rest of the code how to run. So, I'm going to create a case object. I create a case class: </p><precode language="python" precodenum="1"></precode><p>This class, in short, defines case objects through common variables: number of dimensions, inlet and outlet position, and the fluid and porous medium used and their properties. The fluid and porous medium are both represented thorugh dictionaries which other objects or methods can refer to. Let's use it: </p><precode language="python" precodenum="2"></precode><p> So, we've got the same result as the previous, simpler code. We've also created a case object that the other to-be-created code can use. For now, let me compress the code by initializing the superficial velocity in the instantiation (__init__ method):</p><precode language="python" precodenum="3"></precode><p> Later, I'll be able to pass this to a mesh object and Darcy's law method to create a solution that I can plot and output to CSV. </p><span class="meta"><time datetime="2021-04-18T00:00:00+03:00">April 18, 2021</time> · <!--<a href="/tag/css">css</a>--></span></section></main>

