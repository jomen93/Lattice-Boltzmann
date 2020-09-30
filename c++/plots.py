import numpy as np 

x = np.loadtxt("x.dat", unpack=True)
ux = np.loadtxt("ux.dat", unpack=True)

print(np.shape(ux))