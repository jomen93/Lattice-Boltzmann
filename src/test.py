import user
import numpy as np
import matplotlib.pyplot as plt

path = "/Users/johan/Documents/UNAM/Materias/Metodos_numericos/Lattice_Boltzmann/DATA"
ux = np.loadtxt(path+"/ux.dat")*1.22625
uy = np.loadtxt(path+"/uy.dat")*1.22625


def u_theory(nu, dpdy, rho, g, x, d):
    eta = rho*nu
    return (0.5/eta)*(dpdy + rho*g)*((x-0.5*d)**2 - (0.5*d)**2)


nu = user.nu
dpdy = user.dpdy
g = user.g
d = user.d
rho = user.Rho
x = np.linspace(0, d, 100)
print(nu*rho)

plt.plot(x, u_theory(nu, dpdy, rho, g, x, d), "b-", label="theory")
# plt.plot(uy[:, int(len(uy)/2)])
plt.grid(True)
plt.legend()
plt.show()
