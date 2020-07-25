import numpy as np
import matplotlib.pyplot as plt

Nx = 256
Ny = 256

x = np.linspace(0,Nx,Nx)
y = np.linspace(0,Ny,Ny)

x,y = np.meshgrid(x,y)

def sinc(x,y):
	return 0.0001*(x-Nx/2.)**2+0.0001*(y-Ny/2.)**2

plt.imshow(sinc(x,y))
plt.colorbar()
plt.show()
