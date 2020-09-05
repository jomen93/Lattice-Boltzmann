from LB import LatticeBoltzmann
import numpy as np
from utils import Animation
import matplotlib.pyplot as plt

Nx = 64
Ny = 128
tfin = 1000
system = LatticeBoltzmann(Nx, Ny)

system.Init()
ims_vel = []
ims_den = []
for t in range(tfin):
    system.VonKarman_Macroscopic()
    system.Feq_fluids()
    system.Collision()
    system.Streaming()
    system.VonKarman_boundaries()
    if (t % 50 == 0):
        ims_den.append(system.rho)
        ims_vel.append(np.sqrt(system.u[0]**2+system.u[1]**2))
    print("simulation_progess = {:.2f}".format((t/tfin)*100)+"%\r", end="")
Animation(ims_vel, "velocity_karman", plt.cm.viridis)
print("vectorial field done!")
Animation(ims_den, "density_karman", plt.cm.magma)
print("density field done!")
print("Aniamtion_done")
