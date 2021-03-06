# Lattice-Boltzmann

## Navier Stokes equations

<p align = "justify">
A fluid is a continuous medium in space, such a medium is thought of as a composition of small point particles. Classical mechanics is a branch of physics that attempts to describe the behavior of bodies, whether they are solid, liquid or gaseous. It builds its mathematical formalism between experimentation and theory, in fluids it is about building a theory that serves as a model for a subset of real phenomena, due to the nature of the model, it renounces the claim of accuracy in the description and part of the fact. to reflect and allow us to intuit the underlying physical reality.
</p>

<p align="justify">
Fluids in the first approximation are described by the Euler equation, which in summary is to impose the conservation of mass, momentum and energy in the fluid, initially it is built considering that it is a reversible process, since energy losses are not assumed, Euler's equation is written
</p>

<div align = "center">
<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;(\rho&space;v_{i})}{\partial&space;t}&space;=&space;-\frac{\partial&space;\Pi_{ik}}{\partial&space;x_{k}}" title="\frac{\partial (\rho v_{i})}{\partial t} = -\frac{\partial \Pi_{ik}}{\partial x_{k}}" />
</div>

where 

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?\Pi_{ik}&space;=&space;p\delta_{ik}&space;&plus;&space;\rho&space;v_{i}v_{k}" title="\Pi_{ik} = p\delta_{ik} + \rho v_{i}v_{k}" />
</div>

<p align="justify">
This quantity represents a momentum transfer (reversible process) given by the mechanical transport of the fluid particles in space. The starting point to obtain the dynamics of Navier Stokes will be to consider that the viscosity is due to a moment transfer where energy is not conserved, therefore there is an irreversible process, from places where the speed is great to others where speed is small. Then in the definition of the impulse tensor a term is added that accounts for the transfer of viscous impulse. Thus
</p>

<p aling="center">
<img src="https://latex.codecogs.com/gif.latex?\Pi_{ik}&space;=&space;p\delta_{ik}&space;&plus;&space;\rho&space;v_{i}v_{k}&space;-&space;\sigma^{'}_{ik}=-\sigma_{ik}&plus;\rho&space;v_{i}v_{k}" title="\Pi_{ik} = p\delta_{ik} + \rho v_{i}v_{k} - \sigma^{'}_{ik}=-\sigma_{ik}+\rho v_{i}v_{k}" />
</p>

<p align="justify">
defining the stress tensor <img src="https://latex.codecogs.com/gif.latex?\sigma_{ik}" title="\sigma_{ik}" /> and similarly the viscosity stress tensor <img src="https://latex.codecogs.com/gif.latex?\sigma^{'}_{ik}" title="\sigma^{'}_{ik}" />,the latter is considered to express the part of the moment that is not transferred as moment, that is, only by internal friction processes. To find the shape of the viscosity stress tensor, the foundation of fluid particles is studied. The only way that there is friction is to assume that between the fluid particles they have different speeds, therefore there is a relative movement of them, then two approximations are made; in the first approximation, it is assumed that <img src="https://latex.codecogs.com/gif.latex?\sigma^{'}_{ik}" title="\sigma^{'}_{ik}" /> depends on the spatial derivatives of the velocity, in addition it only depends on the first derivatives and its relationship is linear, since <img src="https://latex.codecogs.com/gif.latex?\sigma^{'}_{ik}" title="\sigma^{'}_{ik}" /> must be overridden for <img src="https://latex.codecogs.com/gif.latex?v&space;=&space;\text{cte}" title="v = \text{cte}" />. The second condition that is imposed is that  <img src="https://latex.codecogs.com/gif.latex?\sigma^{'}_{ik}" title="\sigma^{'}_{ik}" /> it must be canceled when there is uniform rotation, since in such a situation there is no internal friction in the fluid. The tensor that satisfies these conditions has the following form
</p>

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?\sigma^{'}_{ik}&space;=&space;\eta\left(\frac{\partial&space;v_{i}}{\partial&space;x_{k}}&space;&plus;&space;\frac{\partial&space;v_{k}}{\partial&space;x_{i}}&space;-&space;\frac{2}{3}\delta_{ik}\frac{\partial&space;v_{l}}{\partial&space;x_{l}}\right)&space;&plus;&space;\zeta\delta_{ik}\frac{\partial&space;v_{l}}{\partial&space;x_{l}}" title="\sigma^{'}_{ik} = \eta\left(\frac{\partial v_{i}}{\partial x_{k}} + \frac{\partial v_{k}}{\partial x_{i}} - \frac{2}{3}\delta_{ik}\frac{\partial v_{l}}{\partial x_{l}}\right) + \zeta\delta_{ik}\frac{\partial v_{l}}{\partial x_{l}}" />
</div>

<p align="justify">
defining the viscous coefficients <img src="https://latex.codecogs.com/gif.latex?\eta&space;=&space;\eta&space;(T,&space;\rho,&space;P)" title="\eta = \eta (T, \rho, P)" /> and <img src="https://latex.codecogs.com/gif.latex?\zeta&space;=&space;\zeta&space;(T,&space;\rho,&space;p)" title="\zeta = \zeta (T, \rho, p)" />, the fluids that depend in a linear way with the strain are called ideal fluids. Then the equation of motion that we have is
</p>

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?\rho\left(\frac{\partial&space;u_{i}}{\partial&space;t}&space;&plus;&space;u_{k}\frac{\partial}{\partial&space;x_{k}}u_{i}\right)&space;=&space;-\frac{\partial&space;p}{\partial&space;x_{i}}&space;&plus;&space;\frac{\partial&space;\sigma^{'}_{ik}}{\partial&space;x_{k}}" title="\rho\left(\frac{\partial u_{i}}{\partial t} + u_{k}\frac{\partial}{\partial x_{k}}u_{i}\right) = -\frac{\partial p}{\partial x_{i}} + \frac{\partial \sigma^{'}_{ik}}{\partial x_{k}}" />
</div>

then

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\sigma^{'}_{ik}}{\partial&space;x_{k}}&space;=&space;\eta\left(\frac{\partial^{2}&space;v_{i}}{\partial&space;x_{k}\partial&space;x_{k}}&plus;&space;\frac{\partial^{2}&space;v_{k}}{\partial&space;x_{k}\partial&space;x_{i}}-\frac{2}{3}\delta_{ik}&space;\frac{\partial^{2}&space;v_{l}}{\partial&space;x_{k}&space;\partial&space;x_{l}}\right)&plus;\zeta&space;\frac{\partial^{2}&space;v_{l}}{\partial&space;x_{k}\partial&space;x_{l}}" title="\frac{\partial \sigma^{'}_{ik}}{\partial x_{k}} = \eta\left(\frac{\partial^{2} v_{i}}{\partial x_{k}\partial x_{k}}+ \frac{\partial^{2} v_{k}}{\partial x_{k}\partial x_{i}}-\frac{2}{3}\delta_{ik} \frac{\partial^{2} v_{l}}{\partial x_{k} \partial x_{l}}\right)+\zeta \frac{\partial^{2} v_{l}}{\partial x_{k}\partial x_{l}}" />
</div>

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\sigma^{'}_{ik}}{\partial&space;x_{k}}=\eta\left(\frac{\partial^{2}&space;v_{i}}{\partial&space;x_{k}\partial&space;x_{k}}&space;&plus;&space;\frac{\partial}{\partial&space;x_{i}}\left(\frac{\partial&space;v_{k}}{\partial&space;x_{k}}\right)-\frac{2}{3}\delta_{ik}\frac{\partial}{\partial&space;x_{k}}\left(\frac{\partial&space;v_{l}}{\partial&space;x_{l}}\right)\right)&space;&plus;&space;\zeta&space;\frac{\partial^{2}&space;v_{l}}{\partial&space;x_{k}\partial&space;x_{l}}" title="\frac{\partial \sigma^{'}_{ik}}{\partial x_{k}}=\eta\left(\frac{\partial^{2} v_{i}}{\partial x_{k}\partial x_{k}} + \frac{\partial}{\partial x_{i}}\left(\frac{\partial v_{k}}{\partial x_{k}}\right)-\frac{2}{3}\delta_{ik}\frac{\partial}{\partial x_{k}}\left(\frac{\partial v_{l}}{\partial x_{l}}\right)\right) + \zeta \frac{\partial^{2} v_{l}}{\partial x_{k}\partial x_{l}}" />
</div>

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\sigma^{'}_{ik}}{\partial&space;x_{k}}=\eta\left(\frac{\partial^{2}&space;v_{i}}{\partial&space;x_{k}\partial&space;x_{k}}\right)&space;&plus;&space;\frac{1}{3}\eta\frac{\partial}{\partial&space;x_{i}}\left(\frac{\partial&space;v_{l}}{\partial&space;x_{l}}\right)&plus;&space;\zeta&space;\frac{\partial^{2}&space;v_{l}}{\partial&space;x_{k}\partial&space;x_{l}}" title="\frac{\partial \sigma^{'}_{ik}}{\partial x_{k}}=\eta\left(\frac{\partial^{2} v_{i}}{\partial x_{k}\partial x_{k}}\right) + \frac{1}{3}\eta\frac{\partial}{\partial x_{i}}\left(\frac{\partial v_{l}}{\partial x_{l}}\right)+ \zeta \frac{\partial^{2} v_{l}}{\partial x_{k}\partial x_{l}}" />
</div>

the equation that describes a fluid in general is given by

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?\rho\left(\frac{\partial&space;u_{i}}{\partial&space;t}&space;&plus;&space;u_{k}\frac{\partial}{\partial&space;x_{k}}u_{i}\right)&space;=&space;-\frac{\partial&space;p}{\partial&space;x_{i}}&space;&plus;&space;\eta\frac{\partial^{2}&space;v_{i}}{\partial&space;x_{k}\partial&space;x_{k}}&space;&plus;&space;\left(\zeta&space;&plus;&space;\frac{1}{3}\eta\right)\frac{\partial}{\partial&space;x_{i}}\frac{\partial&space;v_{l}}{\partial&space;x_{l}}" title="\rho\left(\frac{\partial u_{i}}{\partial t} + u_{k}\frac{\partial}{\partial x_{k}}u_{i}\right) = -\frac{\partial p}{\partial x_{i}} + \eta\frac{\partial^{2} v_{i}}{\partial x_{k}\partial x_{k}} + \left(\zeta + \frac{1}{3}\eta\right)\frac{\partial}{\partial x_{i}}\frac{\partial v_{l}}{\partial x_{l}}" />
</div>

which in a vector form is written as follows

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?\rho&space;\left(\frac{\partial&space;\vec{v}}{\partial&space;t}&plus;(\vec{v}&space;\cdot&space;\nabla&space;)\vec{v}\right)&space;=&space;-\nabla&space;p&space;&plus;\eta&space;\nabla^{2}\vec{v}&space;&plus;&space;\left(\zeta&space;&plus;&space;\frac{1}{3}\eta\right)\nabla(\nabla&space;\cdot&space;\vec{v})" title="\rho \left(\frac{\partial \vec{v}}{\partial t}+(\vec{v} \cdot \nabla )\vec{v}\right) = -\nabla p +\eta \nabla^{2}\vec{v} + \left(\zeta + \frac{1}{3}\eta\right)\nabla(\nabla \cdot \vec{v})" />
</div>

<p align="justify">
when the fluid is considered incompressible, the fact of <img src="https://latex.codecogs.com/gif.latex?\nabla&space;\cdot&space;\vec&space;{v}&space;=&space;0" title="\nabla \cdot \vec {v} = 0" /> is assumed a priori
</p>

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\vec{v}}{\partial&space;t}&plus;(\vec{v}&space;\cdot&space;\nabla&space;)\vec{v}&space;=&space;-\frac{1}{\rho}\nabla&space;p&space;&plus;\frac{\eta}{\rho}&space;\nabla^{2}\vec{v}" title="\frac{\partial \vec{v}}{\partial t}+(\vec{v} \cdot \nabla )\vec{v} = -\frac{1}{\rho}\nabla p +\frac{\eta}{\rho} \nabla^{2}\vec{v}" />
</div>
