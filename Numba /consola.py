import numpy as np 
from numba import jitclass
from numba import int32, float32, float64

spec = [
	("x", int32),
	("y", float32),
	("z", float64[:])]


@jitclass(spec)
class vec(object):
		"""docstring for vec"""
		def __init__(self, x,y):
			self.x = x
			self.y = y
			self.z = np.ones(100)

		def add(self, dx,dy):
			self.x += dx
			self.y += dy


clase = vec(1,2)
print(clase.x) 
print(clase.z)