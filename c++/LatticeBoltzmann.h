const int Lx =7;
const int Ly =7;
const int Q=9;
const int D=2;

class LatticeBoltzmann
{
private:
	double w[Q];
	int cx[Q];
	int cy[Q];
	double f[Lx][Ly][Q];
	double f_post[Lx][Ly][Q];
	double rho[Lx][Ly];
	double ux[Lx][Ly];
	double uy[Lx][Ly];
	double tau;
	double dt;
	
public:
		LatticeBoltzmann(void);

		// Public functions 
        void begin(void);
        void Collision(void);
        void Advection(void);
        void Macroscopic(void);
        void State(void);

        // Statistical distributions
        double feq(int i, double rhoo, double uxo, double uyo);
};
