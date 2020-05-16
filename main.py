from LB import LatticeBoltzmann
import numpy as np 
from utils import plot_3d, Animation
import matplotlib.pyplot as plt

Nx = 64
Ny = 64


x = LatticeBoltzmann(Nx,Ny)

x.Init()
ims = []
#plot_3d(x.rho,"rho")
for i in range(100):
	x.Macroscopic()
	x.Feq()
	x.Collision()
	x.Streaming()
	x.Boundaries()
	ims.append(x.rho)

Animation(ims)