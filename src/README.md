# Lattice Boltzmann

Implementation of method

## Basic Theory

<div align="justify">
	The macroscopic dynamics of the fluid can be viewed as the result of a collective behavior of the microscopic particles in the study system. Naturally, the formalism that adequately describes this will be the Navier Stokes equations with conservation considerations, in order to arrive at the macroscopic equations from the Lattice Boltzmann method, an equilibrium disturbance is used. To begin to introduce this formalism more clearly, the dimensionless Boltzmann equation should be considered in order to better understand the scales on which our system is positioned. We are going to denote the dimensionless quantities with <img src="https://latex.codecogs.com/gif.latex?\widetilde{x}" title="\widetilde{x}" />  and an appropriate characteristic number, we define
</div>

<p align="center">
	<img src="https://latex.codecogs.com/gif.latex?\widetilde{t}&space;=&space;\frac{\xi_{o}}{x_{o}}&space;t&space;\qquad&space;\widetilde{x}&space;=&space;\frac{x}{x_{o}}\qquad&space;\widetilde{\xi}=\frac{\xi}{\xi_{o}}\nonumber\\&space;\widetilde{\tau}&space;=&space;\frac{\xi_{o}}{\overline{x}}\tau\qquad&space;\widetilde{F}&space;=&space;\frac{x_{o}}{\rho\xi_{o}^{2}}F&space;\qquad&space;\widetilde{f}&space;=&space;\frac{c_{o}^{3}}{\rho}f" title="\widetilde{t} = \frac{\xi_{o}}{x_{o}} t \qquad \widetilde{x} = \frac{x}{x_{o}}\qquad \widetilde{\xi}=\frac{\xi}{\xi_{o}}\nonumber\\ \widetilde{\tau} = \frac{\xi_{o}}{\overline{x}}\tau\qquad \widetilde{F} = \frac{x_{o}}{\rho\xi_{o}^{2}}F \qquad \widetilde{f} = \frac{c_{o}^{3}}{\rho}f" />
</p>

<div align="justify">
	When the previous definitions on the Boltzmann equation are used, it is found that in terms of the dimensionless variables, the following dimensionless Boltzmann equation is reached
</div>

<p align="center">
	<img src="https://latex.codecogs.com/gif.latex?\frac{\overline{x}}{x_{o}}\left(\frac{\partial&space;\widetilde{f}}{\partial&space;\widetilde{t}}&plus;\widetilde{\xi}_{\alpha}\frac{\partial&space;\widetilde{f}}{\partial&space;\widetilde{x}_{\alpha}}&plus;\widetilde{F}_{\alpha}\frac{\partial&space;\widetilde{f}}{\partial&space;\widetilde{\xi}_{\alpha}}\right)=-\frac{1}{\widetilde{\tau}}\left(\widetilde{f}-\widetilde{f}^{(0)}\right)" title="\frac{\overline{x}}{x_{o}}\left(\frac{\partial \widetilde{f}}{\partial \widetilde{t}}+\widetilde{\xi}_{\alpha}\frac{\partial \widetilde{f}}{\partial \widetilde{x}_{\alpha}}+\widetilde{F}_{\alpha}\frac{\partial \widetilde{f}}{\partial \widetilde{\xi}_{\alpha}}\right)=-\frac{1}{\widetilde{\tau}}\left(\widetilde{f}-\widetilde{f}^{(0)}\right)" />
</p>


<div align="justify">
	where <img src="https://latex.codecogs.com/gif.latex?\overline{x}" title="\overline{x}" /> is the mean free path and <img src="https://latex.codecogs.com/gif.latex?\overline{t}&space;=&space;\overline{x}/\xi_{o}" title="\overline{t} = \overline{x}/\xi_{o}" /> the mean free time. If <img src="https://latex.codecogs.com/gif.latex?\text{Kn}&space;\rightarrow&space;0" title="\text{Kn} \rightarrow 0" /> both sides of the equation go to zero, the left side by the definition of th Knudsen number, on the other hand, the left side vanishes because the distribution function is very similar to the function balance. We see that as Kn gets bigger the system deviates from equilibrium, then a small paramater <img src="https://latex.codecogs.com/gif.latex?\epsilon" title="\epsilon" /> in introduced to make an expansion and to be able to disturb the Boltzmann equation, the we introduce
</div>



<p align="center">
	<img src="https://latex.codecogs.com/gif.latex?f_{i}&space;=&space;\sum_{n=0}^{&plus;\infty}\epsilon^{n}f_{i}^{(n)}\qquad&space;\frac{\partial}{\partial&space;t}&space;=&space;\sum_{n=0}^{&plus;\infty}\epsilon^{n}\frac{\partial}{\partial&space;t_{n}}\qquad\frac{\partial}{\partial&space;x_{\alpha}}&space;=&space;\sum_{n=0}^{&plus;\infty}\epsilon^{n}\frac{\partial}{\partial&space;x_{\alpha}^{n}}\qquad\frac{\partial}{\partial&space;\xi_{\alpha}}&space;=&space;\sum_{n=0}^{&plus;\infty}\epsilon^{n}\frac{\partial}{\partial&space;\xi_{\alpha}^{n}}" title="f_{i} = \sum_{n=0}^{+\infty}\epsilon^{n}f_{i}^{(n)}\qquad \frac{\partial}{\partial t} = \sum_{n=0}^{+\infty}\epsilon^{n}\frac{\partial}{\partial t_{n}}\qquad\frac{\partial}{\partial x_{\alpha}} = \sum_{n=0}^{+\infty}\epsilon^{n}\frac{\partial}{\partial x_{\alpha}^{n}}\qquad\frac{\partial}{\partial \xi_{\alpha}} = \sum_{n=0}^{+\infty}\epsilon^{n}\frac{\partial}{\partial \xi_{\alpha}^{n}}" />
</p>


<div align="justify">
	the first step to be able to discretize the Boltzmann equation is to be able to put the velocity in discrete vectors, for this we use a projection on the hermite base of the Boltzmann distribution function. We find that the discrete equilibrium function has the following form
</div>


<p align="center">
	<img src="https://latex.codecogs.com/gif.latex?\boxed{&space;f_{i}^{(0)}&space;=&space;\rho&space;w_{i}\left(1&plus;\frac{\xi_{i\alpha}u_{\alpha}}{c_{o}^{2}}&plus;\frac{\xi_{i\alpha}\xi_{i\beta}u_{\alpha}u_{\beta}}{2c_{o}^{4}}-\frac{u_{\alpha}u_{\alpha}}{2c_{o}^{2}}\right)&space;}" title="\boxed{ f_{i}^{(0)} = \rho w_{i}\left(1+\frac{\xi_{i\alpha}u_{\alpha}}{c_{o}^{2}}+\frac{\xi_{i\alpha}\xi_{i\beta}u_{\alpha}u_{\beta}}{2c_{o}^{4}}-\frac{u_{\alpha}u_{\alpha}}{2c_{o}^{2}}\right) }" />
</p>


<div align="justify">
	In this case we consider the classical part of the fluids, because we use the discrete equation without a forcing term, so we are working on a formalism to be able to simulate the Navier-Stokes equations
</div>

<p align="center">
	<img src="https://latex.codecogs.com/gif.latex?\boxed{\frac{\partial&space;f_{i}}{\partial&space;t}&plus;\xi_{i\alpha}\frac{\partial&space;f_{i}}{\partial&space;x_{\alpha}}=-\frac{1}{\tau}\left(f_{i}-f_{i}^{(0)}\right)}" title="\boxed{\frac{\partial f_{i}}{\partial t}+\xi_{i\alpha}\frac{\partial f_{i}}{\partial x_{\alpha}}=-\frac{1}{\tau}\left(f_{i}-f_{i}^{(0)}\right)}" />
</p>

<div align="justify">
	Ahora para poder tener la parte discreta de la anterior ecuación se usa el teorema fundamental del cálculo, obteniendo 
</div>

<p align="center">
	<img src="https://latex.codecogs.com/gif.latex?f_{i}(\vec{x}&plus;\vec{\xi}\Delta&space;t,&space;t&space;&plus;&space;\Delta&space;t)-f_{i}(\vec{x},t)&space;=-&space;\frac{1}{\tau}\int_{0}^{\Delta&space;t}[f_{i}(\vec{x}&plus;\vec{\xi}\lambda,&space;t&space;&plus;&space;\lambda)-f_{i}^{(0)}(\vec{x}&plus;\vec{\xi}\lambda,&space;t&space;&plus;&space;\lambda)]d\lambda" title="f_{i}(\vec{x}+\vec{\xi}\Delta t, t + \Delta t)-f_{i}(\vec{x},t) =- \frac{1}{\tau}\int_{0}^{\Delta t}[f_{i}(\vec{x}+\vec{\xi}\lambda, t + \lambda)-f_{i}^{(0)}(\vec{x}+\vec{\xi}\lambda, t + \lambda)]d\lambda" />
</p>

<div align="justify">
	For the calculations used in these studies, the first order of discretization is used, then the following equation is used
</div>

<p align="center">
	<img src="https://latex.codecogs.com/gif.latex?f_{i}(\vec{x}&plus;\vec{\xi}\Delta&space;t,&space;t&space;&plus;&space;\Delta&space;t)-f_{i}(\vec{x},t)&space;=-&space;\frac{1}{\tau}[f_{i}(\vec{x},&space;t)-f_{i}^{(0)}(\vec{x},&space;t)]" title="f_{i}(\vec{x}+\vec{\xi}\Delta t, t + \Delta t)-f_{i}(\vec{x},t) =- \frac{1}{\tau}[f_{i}(\vec{x}, t)-f_{i}^{(0)}(\vec{x}, t)]" />
</p>

<div align="justify">
	the result is explicit in time, the value of the distribution function in the next step depends on the current step.
</div>





