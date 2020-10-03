// my_class.cpp
#include "LatticeBoltzmann.h" 
#include <iostream> 
#include <fstream>
#include <cstring>
#include<time.h>

using namespace std;

ofstream globalfile;
ofstream history;
char datafile[11] = "global.dat";
char data_history[13] = "history.dat";
time_t ltime;

// Lattice Boltzmann Class implementation
// Basic constructor
LatticeBoltzmann::LatticeBoltzmann(void)
{
     // Set model parameters
     // weights
     w[0]=4.0/9.0;
     w[1] = w[2] = w[3] = w[4] = 1.0/9.0;
     w[5] = w[6] = w[7] = w[8] = 1.0/36.0;

     // grid velocities
     cx[0] = 0; cy[0] = 0;
     cx[1] = 1; cy[1] = 0;
     cx[2] = 0; cy[2] = 1;
     cx[3] =-1; cy[3] = 0;
     cx[4] = 0; cy[4] =-1;
     cx[5] = 1; cy[5] = 1;
     cx[6] =-1; cy[6] = 1;
     cx[7] =-1; cy[7] =-1;
     cx[8] = 1; cy[8] =-1;

     // Parameters
     tau = 0.53;
     dt = 1.0;
     rhoo = 1.0;
     Uxo = 1.0;
     Uyo = 1.0;
}

void LatticeBoltzmann::begin(void)
{
     for (int ix = 0; ix < Lx; ++ix) for (int iy = 0; iy < Ly; ++iy)
     {
               rho[ix][iy] = rhoo;
               ux[ix][iy] = Uxo;
               uy[ix][iy] = Uyo;
               for (int q = 0; q < Q; ++q)
               {
                    f[ix][iy][q] = feq(q, rho[ix][iy], ux[ix][iy], ux[ix][iy]);
               }
     }
}
               
double LatticeBoltzmann::feq(int q, double rhoo, double Uxo, double Uyo)
{
     double cu;
     double U2;
     cu = cx[q]*Uxo + cy[q]*Uyo;
     U2 = Uxo*Uxo + Uyo*Uyo;
     return w[q]*rhoo*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*U2);
}

// double LatticeBoltzmann::rho(int ix, int iy)
// {
//      double suma=0;
//      for (int q = 0; q < Q; ++q)
//           suma += f[ix][iy][q];
//      return suma;
// }


void LatticeBoltzmann::Collision()
{
     double FEQ;
     for (int ix = 0; ix < Lx + 1; ++ix) for (int iy = 0; iy < Ly + 1; ++iy) for (int q = 0; q < Q; ++q)
     {
          FEQ = feq(q, rho[ix][iy],ux[ix][iy], uy[ix][iy]);
          f_post[ix][iy][q] = (1-(dt/tau))*f[ix][iy][q] + (dt/tau)*FEQ;
     }
}

void LatticeBoltzmann::Advection()
{
     int jd;
     int id;
     for (int ix = 0; ix < Lx + 1; ++ix) for (int iy = 0; iy < Ly + 1; ++iy) for (int q = 0; q < Q; ++q)
     {
          jd = iy-cy[q];
          id = ix-cx[q];
          if (jd >= 0 && jd <= Ly && id >= 0 &&id <= Lx)
          {
               f[ix][iy][q] = f_post[jd][id][q];
          }
     }
}

void LatticeBoltzmann::Macroscopic(void)
{
     for (int ix = 0; ix < Lx; ++ix) for (int iy = 0; iy < Ly; ++iy) 
     {
          for (int q = 0; q < Q; ++q)
          {
               rho[ix][iy] = f[ix][iy][q];
          }
          ux[ix][iy] = (f[ix][iy][1] + f[ix][iy][5] + f[ix][iy][8]-
                        f[ix][iy][3] - f[ix][iy][6] - f[ix][iy][7])/rho[ix][iy];
          uy[ix][iy] = (f[ix][iy][5] + f[ix][iy][6] + f[ix][iy][2]-
                        f[ix][iy][7] - f[ix][iy][8] - f[ix][iy][4])/rho[ix][iy];
     }   
}


void LatticeBoltzmann::print_config()
{
     globalfile.open(datafile);
     if (!globalfile.is_open())
     {
          cerr << "Error opening file " << endl;
          exit(EXIT_FAILURE);
     }
     printf("Opened global log file ----> %s\n", datafile);
     ltime = time(NULL);
     globalfile << "************************************************\n";
     globalfile << "*********** Lattice Boltzmann Model ************\n";
     globalfile << "************************************************\n";
     globalfile << "\n";
     globalfile << "Autor: Johan S. Méndez\n";
     globalfile << "Date: " << asctime(localtime(&ltime));
     globalfile << "\n";
     globalfile << "Grid system: D2Q9\n";
     globalfile << "Grid size = " << Lx << " x " << Ly << endl;
     globalfile << "Relaxation time tau = " << tau <<  endl;
     globalfile << "Initial velocity ux = " << Uxo << endl;
     globalfile << "Initial velocity uy = " << Uyo << endl;
     globalfile << "Intial density rho_o = " << rhoo << endl;
     globalfile.close();
}

void LatticeBoltzmann::Report(void)
{
     // Report Progress
     history.open(data_history, fstream::in | fstream::out | fstream::app);
     if (!history.is_open())
     {
          cerr << "Error opening file " << endl;
          exit(EXIT_FAILURE);
     }
     history << "===========================================\n";
     history << "iteration " << it <<" complete !\n";
     history << "elapsed time " << elapsed << endl;
     history << "Physical time \n";
     history << "time step \n";

     history.close();
}

void LatticeBoltzmann::Output(void)
{
}




