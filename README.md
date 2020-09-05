# Lattice-Boltzmann
 
## Navier Stokes equations

A fluid is a continuous medium in space, such a medium is thought of as a composition of small point particles. Classical mechanics is a branch of physics that attempts to describe the behavior of bodies, whether they are solid, liquid or gaseous. It builds its mathematical formalism between experimentation and theory, in fluids it is about building a theory that serves as a model for a subset of real phenomena, due to the nature of the model, it renounces the claim of accuracy in the description and part of the fact. to reflect and allow us to intuit the underlying physical reality.

Fluids in the first approximation are described by the Euler equation, which in summary is to impose the conservation of mass, momentum and energy in the fluid, initially it is built considering that it is a reversible process, since energy losses are not assumed, Euler's equation is written

<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;(\rho&space;v_{i})}{\partial&space;t}&space;=&space;-\frac{\partial&space;\Pi_{ik}}{\partial&space;x_{k}}" title="\frac{\partial (\rho v_{i})}{\partial t} = -\frac{\partial \Pi_{ik}}{\partial x_{k}}" />

where 

<img src="https://latex.codecogs.com/gif.latex?\Pi_{ik}&space;=&space;p\delta_{ik}&space;&plus;&space;\rho&space;v_{i}v_{k}" title="\Pi_{ik} = p\delta_{ik} + \rho v_{i}v_{k}" />

This quantity represents a momentum transfer (reversible process) given by the mechanical transport of the fluid particles in space. The starting point to obtain the dynamics of Navier Stokes will be to consider that the viscosity is due to a moment transfer where energy is not conserved, therefore there is an irreversible process, from places where the speed is great to others where speed is small. Then in the definition of the impulse tensor a term is added that accounts for the transfer of viscous impulse. Thus

<img src="https://latex.codecogs.com/gif.latex?\Pi_{ik}&space;=&space;p\delta_{ik}&space;&plus;&space;\rho&space;v_{i}v_{k}&space;-&space;\sigma^{'}_{ik}=-\sigma_{ik}&plus;\rho&space;v_{i}v_{k}" title="\Pi_{ik} = p\delta_{ik} + \rho v_{i}v_{k} - \sigma^{'}_{ik}=-\sigma_{ik}+\rho v_{i}v_{k}" />

definiendo el tensor de tensiones <img src="https://latex.codecogs.com/gif.latex?\sigma_{ik}" title="\sigma_{ik}" /> y análogamente el tensor de tensiones de la viscocidad $\sigma^{'}_{ik}$, se considera que este último expresa la parte de momento que no se transfiere como momento, es decir, solamente por procesos de fricción interna. Para encontrar la forma del tensor de tensiones de la viscosidad se estudia el fundamento de las partículas de fluido. La única manera de que exista rozamiento es asumir que entre las partículas de fluido tienen velocidades diferentes, por lo tanto existe un movimiento relativo de las mismas, entonces se hacen dos aproximaciones; en la primera aproximación, se asume que $\sigma^{'}_{ik}$ depende de las derivadas espaciales de la velocidad, además que solo depende de las primeras derivadas y su relación es lineal, puesto que $\sigma^{'}_{ik}$ debe anularse para $v = \text{cte}$. La segunda condición que se impone es que $\sigma^{'}_{ik}$ debe anularse cuando se tiene rotación uniforme, dado que en tal situación no se produce rozamiento interno en el fluido. El tensor que cumple tales condiciones tiene la forma siguiente
