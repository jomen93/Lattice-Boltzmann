
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
		self.c2 = 3/5.


		# Parameteres od lattice D2Q9
		self.w = np.array([4.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/36 ,1.0/36,1.0/36,1.0/36.]) 
		self.v = np.array([[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]])
		self.Q = 9

		# initialization macroscopic variables
		self.rho = np.ones((Nx,Ny))
		self.Jx = np.zeros((Nx,Ny))
		self.Jy = np.zeros((Nx,Ny))

		# definition of matrix distribution
		self.feq = np.ones((self.Q,Nx,Ny))
		self.f = np.ones((self.Q,Nx,Ny))
		# Initialization auxiliar varible to streaming process
		self.f_post = np.ones((self.Q,Nx,Ny))
		
	def Feq(self):
			self.feq[0] = self.rho*(1-3*self.c2*(1-4/9.))
			for i in range(1,self.Q):
				self.feq[i] = 3*self.w[i]*(self.c2*self.rho +self.v[i][0]*self.Jx+self.v[i][1]*self.Jy)
			return self.feq

	def Init(self):
		x = np.linspace(-5, 5, self.Nx)
		y = np.linspace(-5, 5, self.Ny)
		x, y = np.meshgrid(x, y)
		#self.rho = np.ones((self.Nx,self.Ny))
		self.rho = np.exp(-((x)**2+(y)**2))
		self.Jx = np.zeros((self.Nx,self.Ny))
		self.Jy = np.zeros((self.Nx,self.Ny))
		self.f = self.Feq()


	# Set up equilibrium distribution 
	

	# Calculate macroscopic variables
	def Macroscopic(self):
		self.rho = self.f[0]+self.f[1]+self.f[2]+self.f[3]+self.f[4]+self.f[5]+self.f[6]+self.f[7]+self.f[8]
		self.Jx = ((self.v[1][0]*self.f[1]+self.v[5][0]*self.f[5]+self.v[8][0]*self.f[8])
			-(self.v[3][0]*self.f[3]+self.v[6][0]*self.f[6]+self.v[7][0]*self.f[7]))/self.rho
		self.Jy = ((self.v[2][1]*self.f[2]+self.v[5][1]*self.f[5]+self.v[6][1]*self.f[6])
			-(self.v[4][1]*self.f[4]+self.v[7][1]*self.f[7]+self.v[8][1]*self.f[8]))/self.rho
		
	def Collision(self):
		self.f_post = self.Wp*self.f + self.W*self.Feq()

	def Streaming(self):
		pass
		



		
	

		

