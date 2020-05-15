import matplotlib.pyplot as plt
from  mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_3d(data,name):
	fig = plt.figure()
	#ax = plt.axes(projection="3d")
	ax = fig.gca(projection='3d')

	x = np.linspace(-5, 5, np.shape(data)[0])
	y = np.linspace(-5, 5, np.shape(data)[1])
	x, y = np.meshgrid(x, y)
	surf = ax.plot_surface(x,y,data,rstride=1,cstride=1,cmap= "YlGn",edgecolor='none',linewidth = 0.2)
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel(name)
	cbar = fig.colorbar(surf, shrink=0.5, aspect=8)
	cbar.set_label(name)
	plt.show()



