import numpy as np 
import matplotlib.pyplot as plt

nx = 32; ny = 32
q = np.random.rand(nx,ny)
print(q)
def delta(q):
	e = 0.5
	Lx = 1
	return np.where(q <=e*Lx, 0, (0.5/(e*Lx))*(1+np.cos((np.pi*(q))/(e*Lx))))




def f(x):
	e = 0.001; Lx = 1
	return (0.5/(e*Lx))*(1+np.cos((np.pi*(x))/(e*Lx)))


x = np.linspace(-10,10,100,dtype = float)
y = np.linspace(-10,10,100,dtype = float)

x,y = np.meshgrid(x,y)

plt.imshow(f(x))
plt.colorbar()
plt.show()
# plt.imshow(q(y))
plt.show()
