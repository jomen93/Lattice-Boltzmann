import numpy as np


class LatticeBoltzmann(object):
    def __init__(self, Nx, Ny):
        super(LatticeBoltzmann, self).__init__()
        # definitios of adimensional parameters
        self.Lx = 1                             # m
        self.Ly = 1                             # m
        self.T = 1                              # s
        self.Ux = 1                             # m/s
        self.Uy = 1                             # m/s
        self.Rho = 1                            # kg/m2 air density
        self.Tau = 0.9                          # Relaxation time
        self.Nu = 0.1                           # Kinematic Viscocity

        # From this point all variables are adimensional
        self.Nx = int(Nx/self.Lx)
        self.Ny = int(Ny/self.Ly)
        self.dt = 1
        self.tau = self.Tau/self.T
        self.W = self.dt/self.tau
        self.Wp = 1-self.W
        self.c2 = 1/3.

        # Parameteres od lattice D2Q9
        self.D = 2
        self.Q = 9
        self.v = np.array([(x, y) for x in [0, -1, 1] for y in [0, -1, 1]])
        self.w = 1/36.*np.ones(self.Q)
        self.w[np.asarray([np.linalg.norm(ci) < 1.1 for ci in self.v])] = 1./9.
        self.w[0] = 4./9.
        # self.w = np.array([4.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,
        #                    1.0/36 ,1.0/36 ,1.0/9,1.0/36,1.0/36.])

        self.l = np.arange(self.Q)[np.asarray([ci[0] < 0 for ci in self.v])]
        self.r = np.arange(self.Q)[np.asarray([ci[0] > 0 for ci in self.v])]
        self.c = np.arange(self.Q)[np.asarray([ci[0] == 0 for ci in self.v])]
        self.up = np.arange(self.Q)[np.asarray([ci[1] > 0 for ci in self.v])]
        self.low = np.arange(self.Q)[np.asarray([ci[1] < 0 for ci in self.v])]

        # fluid Parameters
        # Reynolds number
        # kinematic viscosity
        # Reynolds Number
        # self.Re = (self.Lx*self.Ux)/self.nu
        # Relaxation parameters
        # self.omega = self.dt/self.tau

        self.Re = 100.0
        # characteristic length
        L = (self.Nx*self.Ny)
        # Velocity in lattice units
        uLB = 0.04
        nulb = uLB*L/self.Re
        # Relaxation parameters
        self.omega = 1.0/(3.*(uLB*L/self.Re)+0.5)
        self.nu = nulb

        # initialization macroscopic variables
        self.rho = np.ones((Nx, Ny))
        self.u = np.zeros((self.D, Nx, Ny))
        self.F = np.zeros((self.D, Nx, Ny))

        # definition of distributions
        self.feq = np.zeros((self.Q, Nx, Ny))
        self.f = np.zeros((self.Q, Nx, Ny))
        self.f_Force = np.zeros((self.Q, Nx, Ny))
        self.f_Source = np.zeros((self.Q, Nx, Ny))

        # Initialization auxiliar varible to streaming process
        self.f_post = np.zeros((self.Q, Nx, Ny))

    def Feq(self):
        self.feq[0] = self.rho*(1-3*self.c2*(1-self.w[0]))
        for i in range(1, self.Q):
            term = self.v[i][0]*self.u[0]+self.v[i][1]*self.u[1]
            self.feq[i, :, :] = 3*self.w[i]*(self.c2*self.rho+term)
        return self.feq

    def Feq_fluids(self):
        cu = 3.0*np.dot(self.v, self.u.transpose(1, 0, 2))
        usqr = 3/2.*(self.u[0]**2+self.u[1]**2)
        for i in range(self.Q):
            self.f[i, :, :] = self.rho*self.w[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
        return self.f

    def Init(self):
        self.f = self.Feq_fluids()

    # Calculate macroscopic variables
    def Macroscopic(self):
        self.rho = np.sum(self.f, axis=0)
        #self.u = np.dot(self.v.transpose(),self.f.transpose((1,0,2)))/self.rho

        # force and source quatities
        #self.rho = np.sum(self.f,axis = 0) +
        # 0.5*self.dt*np.sum(self.Source(),axis = 0)
        F_ext = 0.5*self.dt*self.F/np.sum(self.f, axis=0)
        self.u = np.dot(self.v.transpose(), self.f.transpose((1, 0, 2)))/self.rho + F_ext

    def VonKarman_Macroscopic(self):
        self.rho = np.sum(self.f, axis=0)
        self.u = np.dot(self.v.transpose(), self.f.transpose((1, 0, 2)))/self.rho
        # coordinates of the cylinder
        cx = self.Nx/4
        cy = self.Ny/2
        r = self.Ny/9
        uLB = 0.04
        ly = self.Ny-1.0
        vel = np.fromfunction(lambda d, x, y: (1-d)*uLB*(1.0+1e-4*np.sin(y/ly*2*np.pi)), (2, self.Nx, self.Ny))

        self.u[:, 0, :] = vel[:, 0, :]
        self.rho[0, :] = 1./(1.-self.u[0, 0, :])*(np.sum(self.f[self.c, 0, :], axis=0)+2.*np.sum(self.f[self.l, 0, :], axis=0))

    def Source(self):
        cuS = 3.0*np.dot(self.v, self.u.transpose(1, 0, 2))
        usqrS = 3/2.*(self.u[0]**2+self.u[1]**2)

        r1 = 4
        r2 = 4
        x_c1 = int(3*self.Nx/10)
        y_c1 = int(self.Nx/2)
        x_c2 = int(7*self.Nx/10)
        y_c2 = int(self.Nx/2)

        y1, x1 = np.ogrid[-x_c1:self.Nx-x_c1, -y_c1:self.Ny-y_c1]
        y2, x2 = np.ogrid[-x_c2:self.Nx-x_c2, -y_c1:self.Ny-y_c2]

        mask1 = x1**2 + y1**2 < r1**2
        mask2 = x2**2 + y2**2 < r2**2

        for i in range(self.Q):
            qo1 = 1e-2
            # qo1 = 1*np.cos(t/0.03)**2
            qo2 = 5e-2
            #qo2 = 1*np.cos(t/0.03)**2

            aux = self.w[i]*(1.+cuS[i]+0.5*cuS[i]**2-usqrS)

            self.f_Source[i, mask1] = aux[mask1]*qo1
            self.f_Source[i, mask2] = aux[mask2]*qo2

        return self.f_Source

    def Pouseuille_Force(self):
        g = 1e-6
        self.F[1] = -g*self.rho
        cuF = 3.0*np.dot(self.v, self.F.transpose(1, 0, 2))
        usqrF = 3/2.*(self.F[0]*self.u[0]+self.F[1]*self.u[1])
        for i in range(self.Q):
            self.f_Force[i, :, :] = self.w[i]*(cuF[i]+0.5*cuF[i]**2-usqrF)
        return self.f_Force

    def Collision(self):
        self.f_post = (1-self.omega)*self.f+self.omega*self.Feq_fluids()
        +(1.0-0.5*(self.omega))*self.Pouseuille_Force()
      # +(1.0-0.5*(self.omega))*self.Source()

    def Streaming(self):
        for i in range(self.Q):
            self.f[i, :, :] = np.roll(np.roll(self.f_post[i, :, :],
                                      self.v[i, 0], axis=0),
                                      self.v[i, 1], axis=1)

    #### implementation Boundary conditions in Fluids ###
    # We focused on analitycal and standard numerical procedures
    # for boundary conditions in fluid flow problems.
    #Bounce Back boundary condition

    def Bounce_Back(self):
    # this module changes acord to specific problem and obstacle
        # l wall
        self.f[self.r, :, 0] = self.f_post[self.l, :, 0]
        # r Wall
        self.f[self.l, :, -1] = self.f_post[self.r, :, -1]
        # up Wall
        self.f[self.low, -1, :] = self.f_post[self.up, -1, :]
        # low Wall
        self.f[self.up, 0, :] = self.f_post[self.low, 0, :]

    def Boundaries(self):
    # Free boundary condition
        # r wall
        self.f[self.l, -1, :] = self.f[self.l, -2, :]
        # l wall
        self.f[self.r, 0, :] = self.f[self.r, 1, :]
        # up wall
        self.f[self.low, :, -1] = self.f[self.low, :, -2]
        # low wall
        self.f[self.up, :, 0] = self.f[self.up, :, 1]
        # corners
        # r uppper
        self.f[1, -1, -1] = self.f[1, -2, -2]
        self.f[3, -1, -1] = self.f[3, -2, -2]
        self.f[4, -1, -1] = self.f[4, -2, -2]
        self.f[0, -1, -1] = self.f[0, -2, -2]
        # l up
        self.f[6, 0, -1] = self.f[6, 1, -2]
        self.f[1, 0, -1] = self.f[1, 1, -2]
        self.f[7, 0, -1] = self.f[7, 1, -2]
        self.f[0, 0, -1] = self.f[0, 1, -2]
        # r low
        self.f[2, -1, 0] = self.f[2, -1, 1]
        self.f[3, -1, 0] = self.f[3, -1, 1]
        self.f[5, -1, 0] = self.f[5, -1, 1]
        self.f[0, -1, 0] = self.f[0, -1, 1]
        # l low
        self.f[2, 0, 0] = self.f[2, 1, 1]
        self.f[6, 0, 0] = self.f[6, 1, 1]
        self.f[8, 0, 0] = self.f[8, 1, 1]
        self.f[0, 0, 0] = self.f[0, 1, 1]

    def Pouseuille_Boundaries(self):
        # Bounce Back in top and bottom boundaries to replicate the not slip
        # boundary conditions

        # self.f[self.up,0,:] = self.f_post[self.low,0,:]
        # self.f[self.low,-1,:] = self.f_post[self.up,-1,:]
        # self.f[self.low,0,:] = self.f_post[self.up,0,:]
        # self.f[self.up,-1,:] = self.f_post[self.low,-1,:]

        self.f[self.up, -1, :] = self.f_post[self.low, -1, :]
        self.f[self.low, 0, :] = self.f_post[self.up, 0, :]

        # Periodic boundaries in l and r boundaries
        self.f[self.r, :, 0] = self.f_post[self.r, :, -1]
        self.f[self.l, :, -1] = self.f_post[self.l, :, 0]

        # self.f[6,:,0] = self.f_post[6,:,-1]
        # # corners
        # # r up
        # self.f[8,-1,-1] = self.f_post[4,-1,-1]
        # self.f[6,-1,-1] = self.f_post[3,-1,-1]
        # self.f[2,-1,-1] = self.f_post[1,-1,-1]
        # # l up
        # self.f[2,-1,0] = self.f_post[1,-1,0]
        # self.f[5,-1,0] = self.f_post[7,-1,0]
        # self.f[3,-1,0] = self.f_post[6,-1,0]
        # # r low
        # self.f[7,0,-1] = self.f_post[5,0,-1]
        # self.f[1,0,-1] = self.f_post[2,0,-1]
        # self.f[6,0,-1] = self.f_post[3,0,-1]
        # # l low
        # self.f[3,0,0] = self.f_post[6,0,0]
        # self.f[4,0,0] = self.f_post[8,0,0]
        # self.f[1,0,0] = self.f_post[2,0,0]
        # obstacle
        # Definition of obstacle
        # r = 16
        # x_c1 = int(3*self.Nx/6); y_c1 = int(self.Nx/2)
        # y,x = np.ogrid[-x_c1:self.Nx-x_c1, -y_c1:self.Ny-y_c1]
        # mask = x**2 + y**2 < r**2
        # self.f[:,mask] = self.f_post[:,mask]

    def VonKarman_boundaries(self):
        # r Wall: outflow condition
        self.f[self.r, -1, :] = self.f[self.r, -2, :]
        # l Wall
        self.f[self.r, 0, :] = self.f[self.l, 0, :] + self.Feq_fluids()[self.r, 0, :] - self.f[self.l, 0, :]
        cx = self.Nx/4; cy=self.Ny/2; r=self.Ny/9 # coordinates of the cylinder 
        obstacle = np.fromfunction(lambda x, y: (x-cx)**2+(y-cy)**2 < r**2, (self.Nx, self.Ny))
        noslip = [self.v.tolist().index((-self.v[i]).tolist()) for i in range(self.Q)]
        for i in range(self.Q): 
            self.f_post[i, obstacle] = self.f[noslip[i], obstacle]
