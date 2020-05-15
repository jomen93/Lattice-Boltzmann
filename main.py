from LB import LatticeBoltzmann
import numpy as np 
from utils import plot_3d
Nx = 4
Ny = 4

x = LatticeBoltzmann(Nx,Ny)


x.Init()
x.Macroscopic()
#plot_3d(x.rho,"rho")
x.Collision()
print(x.f_post)

#x.Macroscopic()
#print(x.Colision())

