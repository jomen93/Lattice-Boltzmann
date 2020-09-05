import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy import asarray as ar, exp

ux = np.loadtxt("ux.dat")
uy = np.loadtxt("uy.dat")
rho = np.loadtxt("rho.dat")

u = np.sqrt(ux**2 + uy**2)
lenu = np.shape(u)[0]

x = np.linspace(0,128,100)

print(np.mean(rho))

def uteo(x):
	xbot = 64
	xtop = 0 
	Rho = np.mean(rho)
	nu = 0.003333333333333336
	# return (h**2 - (x-64)**2)
	return -(1/(2*Rho*nu))*((x-xbot)*(x-xtop))



# interpolation of simulation points 
yy = u[:,int(lenu/2)]
xx = np.linspace(0,128,len(yy))

def parabola(x,a,b,c):
    return a*(x-b)*(x-c)

popt,pcov = curve_fit(parabola,xx,yy)

plt.plot(xx,yy,'b.',label='Simulación')
plt.plot(x,parabola(x,*popt),'r-',label='Ajuste')
plt.grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.5)
plt.title("Perfil de Velocidad")
plt.ylabel("$|v|$")
plt.xlabel("Posicion")
plt.savefig("vel_profile",dpi = 1200, transparent = True)
plt.legend()
plt.show()

print(popt)
print(pcov)



# plt.plot(u[:,int(lenu/2)],"b.")
# # plt.plot(x,uteo(x),"r-")
# plt.grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.5)
# plt.title("Perfil de Velocidad")
# plt.ylabel("$|v|$")
# plt.xlabel("Posicion")
# plt.savefig("vel_profile",dpi = 1200)
# plt.show()

