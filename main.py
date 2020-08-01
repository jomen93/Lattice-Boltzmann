from LB import LatticeBoltzmann
import numpy as np 
from utils import plot_3d, Animation
import matplotlib.pyplot as plt

Nx = 512
Ny = 512
tfin = 500

x = LatticeBoltzmann(Nx,Ny)

x.Init()
ims_vel = []
ims_den = []
#plot_3d(x.rho,"rho")
print("Reynolds Number :", x.Re)
for t in range(tfin):
	x.Macroscopic()
	x.Feq_fluids()
	# x.Feq()
	x.Collision()
	x.Streaming()
	# x.Bounce_Back()
	x.Boundaries()
	if (t%5 == 0):
		ims_den.append(x.rho)
		ims_vel.append(np.sqrt(x.J[0]**2+x.J[1]**2))
	print("simulation_progess = {:.2f}".format((t/tfin)*100)+"%\r",end ="")
print("")
print("Simulation_done")
Animation(ims_vel,"velocity_fluids",plt.cm.seismic)
Animation(ims_den,"density_fluids",plt.cm.autumn)
print("Aniamtion_done")