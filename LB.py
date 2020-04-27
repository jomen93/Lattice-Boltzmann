
import numpy as np 

class LatticeBoltzmann(object):
	"""docstring for LatticeBoltzmann"""
	def __init__(self,c,tau,Lx,Ly,model):
		super(LatticeBoltzmann, self).__init__()
		self.c = c
		self.tau = tau
		self.Lx = Lx
		self.Ly = Ly
		self.model =model 

		if model == True:

			self.w = np.array([4.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/36 ,1.0/36,1.0/36,1.0/36.]) 
			self.v = np.array([[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]])
			self.f = np.zeros(9)	
	
	def rho(self):
		return np.sum(self.f)

	def J(self):
		return np.dot(self.f,self.v)


	def feq(self,rho,J):
		"""
		Compute de equilibrium distribution for the waves	
		"""	
		Feq = np.zeros(9)
		for i in range(len(Feq)):
			if i == 0:
				Feq[i] = rho*(1-(3*self.c**2)*(1-self.w[i]))
			else:
				Feq[i] = 3*self.w[i]*(self.c**2*rho+(np.dot(self.v[i],J)))
		return Feq

	def Initialize(self,rho,J):
		self.f = self.feq(rho,J)
		return self.f

	def Colsion(self):
		pass

	def Streaming(self):
		pass

		

