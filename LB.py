
import numpy as np 

class LatticeBoltzmann(object):
	"""docstring for LatticeBoltzmann"""
	def __init__(self,Nx,Ny):
		super(LatticeBoltzmann, self).__init__()
		self.Nx = Nx
		self.Ny = Ny
		self.dt = 1
		self.tau = 0.5 
		self.W = self.dt/self.tau; 
		self.Wp = 1-self.W
		self.c2 = 2/5.


		# Parameteres od lattice D2Q9
		self.D = 2; self.Q = 9 
		self.v = np.array([(x,y) for x in [0,-1,1] for y in [0,-1,1]])
		self.w = np.array([4.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/36 ,1.0/36 ,1.0/9,1.0/36,1.0/36.]) 

		self.left = np.arange(self.Q)[np.asarray([ci[0]<0  for ci in self.v])] 
		self.center = np.arange(self.Q)[np.asarray([ci[0]==0  for ci in self.v])] 
		self.right = np.arange(self.Q)[np.asarray([ci[0]>0  for ci in self.v])] 
		self.upper = np.arange(self.Q)[np.asarray([ci[1]>0  for ci in self.v])] 
		self.lower = np.arange(self.Q)[np.asarray([ci[1]<0  for ci in self.v])] 
		
		# initialization macroscopic variables
		self.rho = np.ones((Nx,Ny))
		self.J = np.zeros((self.D,Nx,Ny))
		#self.Jy = np.zeros((Nx,Ny))

		# definition of matrix distribution
		self.feq = np.zeros((self.Q,Nx,Ny))
		self.f = np.zeros((self.Q,Nx,Ny))
		# Initialization auxiliar varible to streaming process
		self.f_post = np.zeros((self.Q,Nx,Ny))
		
	def Feq(self):
			self.feq[0] = self.rho*(1-3*self.c2*(1-self.w[0]))
			for i in range(1,self.Q):
				self.feq[i,:,:] = 3*self.w[i]*(self.c2*self.rho +self.v[i][0]*self.J[0]+self.v[i][1]*self.J[1])
			return self.feq

	def Init(self):
		x = np.linspace(0, 5, self.Nx)
		y = np.linspace(0, 5, self.Ny)
		x, y = np.meshgrid(x, y)
		self.rho = np.exp(-10*((x-5/2)**2+(y-5/2)**2))
		self.f = self.Feq()	

	# Calculate macroscopic variables
	def Macroscopic(self):
		self.rho = np.sum(self.f,axis = 0)
		self.J = np.dot(self.v.transpose(),self.f.transpose((1,0,2)))/self.rho
		
	def Collision(self):
		self.f_post = self.Wp*self.f + self.W*self.Feq()

	def Streaming(self):
		for i in range(self.Q):
			self.f[1,:,:] = np.roll(np.roll(self.f_post[i,:,:],self.v[i,0],axis = 0),self.v[i,1],axis=1)

	def Boundaries(self):
		pass
		#self.f[self.right,1:-1,1] = self.f_post[self.right,1:-1,-1] 
		#self.f[self.left,:,-1] = self.f[self.right,:,1] 
		#self.f[self.left,:,-1] = self.f[self.left,:,0]
		
		#self.f[self.upper,-1,:] = self.f[self.upper,0,:]
		#self.f[self.left,:,-1] = self.f[self.left,:,0]

		#self.f[self.right,0,:] = self.f[self.right,1,:]
		#self.f[self.upper,:,-1] = self.f[self.upper,:,-2]
		# self.f[self.upper,:,0] = self.f[self.upper,:,1]
		
		



		
	

		

