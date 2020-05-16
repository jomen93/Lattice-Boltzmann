from LB import LatticeBoltzmann
import numpy as np 
from utils import plot_3d
Nx = 32
Ny = 32


x = LatticeBoltzmann(Nx,Ny)

x.Init()

plot_3d(x.rho,"rho")
for _ in range(5):
	x.Macroscopic()
	x.Feq()
	x.Collision()
	x.Streaming()
	x.Boundaries()
	#print(x.J[1])
	plot_3d(x.rho,"rho")
	

