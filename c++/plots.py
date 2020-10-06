import numpy as np 
import matplotlib.pyplot as plt

ux = np.loadtxt("ux.dat", unpack=True)
uy = np.loadtxt("uy.dat", unpack=True)
rho = np.loadtxt("rho.dat", unpack=True)

plt.imshow(ux)
plt.colorbar()
plt.show()
plt.imshow(uy)
plt.colorbar()
plt.show()
plt.imshow(rho)
plt.colorbar()
plt.show()
