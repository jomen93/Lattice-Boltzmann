#include <iostream> 
#include <fstream> 
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#include "my_class.h"


int main (int argc, char *argv[])
{
	LatticeBoltzmann LB;
	int tmax = 1;
	LB.begin();
	for (int t = 0; t < tmax; ++t)
	{
		LB.Collision();
		LB.Advection();
		LB.Macroscopic();
	}
	LB.State();
	return 0; 
}