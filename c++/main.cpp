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
	LatticeBoltzmann LB;
	time_t start, end; 
	int tmax = 10;
	time(&start);
	double elapsed;
	LB.begin();
	for (int t = 0; t < tmax; ++t)
	{
		LB.Collision();
		LB.Advection();
		LB.Macroscopic();
	}
	LB.State();
	time(&end);
	elapsed = end - start;
	printf("\n **** All runs completed **** \n");
	printf("**** elapesed time = %f s ****", elapsed);
	printf(" Using \n");
	return 0; 
}