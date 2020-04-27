from LB import LatticeBoltzmann
import numpy as np 

wo = 1/3.
c = 0.5
tau = 0.5
Lx = 128
Ly = 128


x = LatticeBoltzmann(c,tau,Lx,Ly,True)
j = np.zeros(2)
rho = 1

print(x.Initialize(rho,j))
print(x.rho())
print(x.J())
