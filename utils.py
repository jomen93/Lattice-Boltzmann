import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
from  mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numpy as np

def plot_3d(data,name):

	X = np.linspace(0, 5, np.shape(data)[0])
	Y = np.linspace(0, 5, np.shape(data)[1])
	x, y = np.meshgrid(X, Y)

	
	mappable = plt.cm.ScalarMappable(cmap = plt.cm.viridis)
	mappable.set_array(data)
	#mappable.set_clim(0.5,1.1)

	fig = plt.figure(figsize = (10,4))
	
	ax1 = fig.add_subplot(121,projection="3d")
	ax1.plot_surface(x,y,data,cmap=mappable.cmap, norm = mappable.norm, linewidth=0,antialiased=False)
	ax1.set_xlabel("x")
	ax1.set_ylabel("y")
	ax1.set_zlabel(name)
	ax1.set_xlim(np.min(X),np.max(X)) 
	ax1.set_ylim(np.min(Y),np.max(Y)) 
	#ax1.set_zlim(0.5,1.1) 
	
	
	ax2 = fig.add_subplot(122)
	ax2.imshow(data,cmap=mappable.cmap, norm=mappable.norm,
		extent=(np.min(x),np.max(x),np.min(y),np.max(y)),interpolation = "none")

	plt.colorbar(mappable)
	plt.tight_layout()
	plt.show()

def Animation(cube):
	fig = plt.figure()
	#Writer = animation.writers['ffmpeg']
	#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

	mappable = plt.cm.ScalarMappable(cmap = plt.cm.jet)
	mappable.set_array(cube[-1])
	mappable.set_clim(-1,1)
	ims = []
	for i in range(len(cube)):
		im = plt.imshow(cube[i],cmap=mappable.cmap, 
			interpolation = "sinc",animated=True)
		ims.append([im])
	ani = animation.ArtistAnimation(fig,ims, interval = 50,blit=True,repeat_delay=1000)
	plt.colorbar(mappable)
	plt.tight_layout()
	#ani.save('im.mp4', writer=writer)
	plt.show()


