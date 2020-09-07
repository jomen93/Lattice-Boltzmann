from LB import LatticeBoltzmann as LB
import numpy as np
from utils import Animation, plot_3d
import matplotlib.pyplot as plt
import user

Nx = user.Nx
Ny = user.Ny
tfin = user.tfin
Rho = user.Rho
Tau = user.Tau

x = LB(Nx, Ny, Rho, Tau)

x.Init()
ims_vel = []
ims_den = []
#plot_3d(x.rho,"rho")
print("Reynolds Number :", x.Re)
for t in range(tfin):
    x.Macroscopic()
    x.Feq_fluids()
    x.Collision()
    x.Streaming()
    x.Pouseuille_Boundaries()
    # Save images each 50 steps
    if (t % 50 == 0):
        ims_den.append(x.rho)
        ims_vel.append(np.sqrt(x.u[0]**2+x.u[1]**2))
    print("simulation_progess = {:.2f}".format((t/tfin)*100)+"%\r", end="")
np.savetxt('ux.dat', x.u[0])
np.savetxt('uy.dat', x.u[1])
np.savetxt('rho.dat', x.rho)
plot_3d(np.sqrt(x.u[0]**2+x.u[1]**2), "norm_velocity", plt.cm.viridis)
print("viscocity in lattice units :", x.Nu)
print("density in lattice units :", np.mean(x.rho))
print("")
print("Simulation_done")
Animation(ims_vel, "velocity_fluids", plt.cm.viridis)
print("vectorial field done!")
Animation(ims_den, "density_fluids", plt.cm.magma)
print("density field done!")
print("Aniamtion_done")
