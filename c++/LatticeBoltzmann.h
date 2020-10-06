#ifndef LATTICE_BOLTZMANN_H
#define LATTICE_BOLTZMANN_H

//  Lattice Boltzmann Class
class LatticeBoltzmann
{
public:
	//  
	int it;
	float elapsed;
	// Lattice Parameters
	static const int Lx =64;
	static const int Ly =64;
	static const int D=2;
	static const int Q=9;
	
	// Lattice Variables 
	double w[Q];
	int cx[Q];
	int cy[Q];
	double f[Lx][Ly][Q];
	double f_post[Lx][Ly][Q];
	double rho[Lx][Ly];
	double ux[Lx][Ly];
	double uy[Lx][Ly];

	// physical parameters
	double tau;
	double dt;
	double rhoo;
	double Uxo;
	double Uyo;
	// Member functions
	LatticeBoltzmann(void);

	// Public functions 
    void begin(void);
    void Collision(void);
    void Advection(void);
    void Macroscopic(void);
    void print_config(void);
    void Report(void);
    void Output(void);

    // Statistical distributions
    double feq(int i, double rhoo, double uxo, double uyo);
};

#endif // LATTICE_BOLTZMANN_H