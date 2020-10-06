#include <iostream> 
#include <fstream> 
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

#include "LatticeBoltzmann.h"


int main (int argc, char *argv[])
{
	time_t start, end;
	int tmax = 2;

	LatticeBoltzmann LB;
	LB.print_config();
	LB.begin();
	for (int t = 0; t < tmax; ++t)
	{
	    start = time(NULL); 
		LB.it = t;
		// Perform the colision step
		LB.Collision();
		// Perform the advection in the grid
		LB.Advection();
		// Calculate of Macroscopic Physical quantities 
		LB.Macroscopic();
		// Print in Log file the parameters of simulation 
		end = time(NULL);
		LB.elapsed = double(end-start);
		LB.Report();
	}
	LB.Output();
	return 0; 
}