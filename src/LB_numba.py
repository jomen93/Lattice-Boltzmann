import numpy as np
from numba import jitclass
from numba import int32, int64, float32, float64
from numba import typeof

spec = [
	("Nx", int32),
	("Ny", int64),
	("dt", int32),
	("tau", float32),
	("W", float32),
	("Wp", float32),
	("c2", float32),
	("D", int32),
	("Q", int32),
	("v", int64[:,:]),
	("w", float64[:]),
	("left", int64[:]),
	# ("right", int32[:]),
	# ("center", int32[:]),
	# ("upper", int32[:]),
	# ("lower", int32[:]),
]

@jitclass(spec)
class LatticeBoltzmann(object):
	def __init__(self,Nx,Ny):
		self.Nx = Nx
		self.Ny = Ny
		self.dt = 1
		self.tau = 0.5 
		self.W = self.dt/self.tau; 
		self.Wp = 1-self.W
		self.c2 = 3/5.

		# # Parameteres od lattice D2Q9
		self.D = 2; self.Q = 9 
		self.v = np.array([(x,y) for x in [0,-1,1] for y in [0,-1,1]])
		self.w = np.array([4.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/36 ,1.0/36 ,1.0/9,1.0/36,1.0/36.]) 

		# self.left = np.arange(self.Q)[np.asarray([ci[0]<0  for ci in self.v])] 
		# self.right = np.arange(self.Q)[np.asarray([ci[0]>0  for ci in self.v])] 
		# self.center = np.arange(self.Q)[np.asarray([ci[0]==0  for ci in self.v])] 
		# self.upper = np.arange(self.Q)[np.asarray([ci[1]>0  for ci in self.v])] 
		# self.lower = np.arange(self.Q)[np.asarray([ci[1]<0  for ci in self.v])]

x = LatticeBoltzmann(5,5)
