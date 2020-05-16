from LB import LatticeBoltzmann
import numpy as np 
from utils import plot_3d
Nx = 128
Ny = 128


x = LatticeBoltzmann(Nx,Ny)

x.Init()


plot_3d(x.rho,"rho")
for _ in range(100):
	x.Macroscopic()
	x.Feq_fluid()
	x.Collision()
	x.Streaming()
	x.Boundaries()

	#print(x.J[1])
	plot_3d(x.rho,"rho")
	

