from LB import LatticeBoltzmann
import numpy as np 
from utils import plot_3d
import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

Nx = 256
Ny = 256


x = LatticeBoltzmann(Nx,Ny)

x.Init()
fig = plt.figure()
ims = []
#plot_3d(x.rho,"rho")
for _ in range(100):
	x.Macroscopic()
	x.Feq_fluid()
	x.Collision()
	x.Streaming()
	x.Boundaries()
	ims.append([plt.imshow(x.rho,animated=True,cmap=plt.cm.jet,interpolation="spline36")])
ani = animation.ArtistAnimation(fig,ims, interval = 50,blit=True,repeat_delay=1000)
#im_ani.save('im.mp4', writer=writer)
plt.show()	

