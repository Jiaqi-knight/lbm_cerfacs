#include <iostream>
#include <time.h>  //clock_t, clock()

#include "./CommonFunctions.hpp"
#include "./Output.hpp"
#include "./InputFile.hpp"
#include "./Memory.hpp"
#include "./ParamFile.hpp"
#include "GlobalVars.hpp"
#include <omp.h>

using namespace std;

// Global variables ################################################
string LBMmodel_("D2Q9");   // choice of LBMmodel_ : D2Q4, D2Q5 or D2Q9 
int D_(2);                  // Dimension
int Q_(0);                  // lattice number of speeds
double *ex_(0), *ey_(0);    // lattice speed coordinates
double *w_(0);              // lattice weights
double cs_(0);              // lattice sound speed
double lx_(0), ly_(0);      // Size of the computational domain in meters
double dx_(0);              // Mesh size
double dt_(0);              // time step
int Nx_(0), Ny_(0);         // Number of cells
int nite_(0);               // Number of time steps (nite_=0: conv criteria)
double simTime_(0.);        // Simulation time (in seconds)
double Nu_(0);              // Kinematic viscosity in physical units
double NuLb_(0);            // Kinematic viscosity in lattice units
double tau_(0);             // Relaxation time (for explicit LBE)
double Tref_(0);            // Reference temperature
double Pref_(0);            // Reference pressure
int freqOut_(1);            // Output frequency
int freqIntQuant_(1000);    // frequency of integrated quantities
int aeroOrder_(2);          // Development order of fEq
double Rgas_(287.15);       // Ideal gas constant (air)
double Gamma_(1.4);         // heat capacity ratio
double cryoRatio_(1.);      // Cryogenic ratio = Mach cryo / real Mach

int nThreads_(1);            //number of threads openMP

double *Hxx_(0), *Hxy_(0), *Hyy_(0);    // Hermite polynomials

bool reg_(false);
string config_file;

double tCPU_(0.);          // Elapsed time (Total)
double tEq_(0.);           // Elapsed time (Equilibrium)
double tReg_(0.);          // Elapsed time (Regularization)
double tColl_(0.);         // Elapsed time (Collision)
double tStream_(0.);       // Elapsed time (Streaming)
double tMacro_(0.);        // Elapsed time (Macroscopic)
double tOut_(0.);          // Elapsed time (Output)
// End of Global variables ###########################################


int main(int argc, char* argv[]){


  // ================ Gestion of argument (input vtk file) =============== //

  int endSimul(0) ;
  int iter(0);                  // Iteration number
  ifstream param_file;
  string config_file_name("param.txt");

  printHeader();

  if (argc == 2) config_file_name =  argv[1]; 
  param_file.open(&config_file_name[0], ifstream::in);
  if (!param_file){
    cout << "\n\tThe configuration file "<< config_file_name << " is missing!\n" << endl;
    exit (EXIT_FAILURE);
  }   

  // ======== Read parameters and initialisation (physical units) ======== //

  readParamFile(param_file); // recover Q_

  getMeshData(config_file);  // recover Nx_, Ny_ and dx_ values

  omp_set_num_threads(nThreads_);


  double **rho = allocateDbl(Nx_, Ny_);
  double **ux  = allocateDbl(Nx_, Ny_);
  double **uy  = allocateDbl(Nx_, Ny_);
  double **ux0 = allocateDbl(Nx_, Ny_);
  double **uy0 = allocateDbl(Nx_, Ny_);

  double ***f       = allocateDbl(Nx_, Ny_, Q_);   
  double ***fEq     = allocateDbl(Nx_, Ny_, Q_);
  double ***fStream = allocateDbl(Nx_, Ny_, Q_);

  initFromConfigFile(config_file, rho, ux, ux0, uy, uy0);

  generateLattice();   

  convertToLattice(rho, ux, uy, ux0, uy0);   // convert to lattice units

  writeOutputVTK(0, rho, ux, uy);

  ofstream IntQuant_file;
  IntQuant_file.open("output/IntegratedQuantities.txt");
  computeIntegratedQuantities(IntQuant_file, rho, ux, uy, 0);

  EquilibriumDF(rho, ux, uy, f);

  // ============================= Time Loop ============================= //
  // set format manipulators 
  cout << fixed << showpoint; 
  
  cout << endl << "\t========== Starting iterations ==========" << endl ;
  clock_t tmp_tCPU = clock();

  while (endSimul==0){

    iter+=1;

    if ((iter%freqOut_) == 0 or iter==nite_){ 
      cout << "\nit = " << right << setw(7) << iter << left << " || time = " << setw(8) << iter * dt_ << " (s)" ;
    }

    EquilibriumDF(rho, ux, uy, fEq);

    if (reg_) regularization(f, fEq); // To be modified

    collisionBGK(f, fEq);

    streaming(f, fStream);

    macroVar(rho, ux, uy, f);

    if ( (iter%freqOut_) == 0 or iter==nite_) writeOutputVTK(iter, rho, ux, uy);

    if ( (iter%freqIntQuant_) == 0 or iter==nite_) computeIntegratedQuantities(IntQuant_file, rho, ux, uy, iter*dt_);

    // End loop : convergence criteria or t=Nite
    if (nite_>0 && iter==nite_ ) endSimul=1;
  }
  cout << "\n\n" << "\t=========== Ending iterations ===========" << "\n\n";
  tmp_tCPU = clock() - tmp_tCPU;
  tCPU_ += (double) tmp_tCPU / CLOCKS_PER_SEC;

  cout<< left << setw(35) <<"\t Total CPU time: "                 << setw(8) << tCPU_    << " (s)   " << setprecision(2) << setw(6) << right << tCPU_    / tCPU_ * 100. << " %\n\n" << setprecision(6);
  cout<< left << setw(35) <<"\t CPU time for Eq operator: "       << setw(8) << tEq_     << " (s)   " << setprecision(2) << setw(6) << right << tEq_     / tCPU_ * 100. << " %\n" << setprecision(6);
  cout<< left << setw(35) <<"\t CPU time for Reg operator: "      << setw(8) << tReg_    << " (s)   " << setprecision(2) << setw(6) << right << tReg_    / tCPU_ * 100. << " %\n" << setprecision(6);
  cout<< left << setw(35) <<"\t CPU time for Stream operator: "   << setw(8) << tStream_ << " (s)   " << setprecision(2) << setw(6) << right << tStream_ / tCPU_ * 100. << " %\n" << setprecision(6);
  cout<< left << setw(35) <<"\t CPU time for Macro operator: "    << setw(8) << tMacro_  << " (s)   " << setprecision(2) << setw(6) << right << tMacro_  / tCPU_ * 100. << " %\n" << setprecision(6);
  cout<< left << setw(35) <<"\t CPU time for Coll operator: "     << setw(8) << tColl_   << " (s)   " << setprecision(2) << setw(6) << right << tColl_   / tCPU_ * 100. << " %\n" << setprecision(6);
  cout<< left << setw(35) <<"\t CPU time for Output operator: "   << setw(8) << tOut_    << " (s)   " << setprecision(2) << setw(6) << right << tOut_    / tCPU_ * 100. << " %\n" << setprecision(6);
  cout<< setprecision(2);
  cout<< left << setw(35) <<"\t CPU time/iteration/point: "<< setw(8) << scientific << tCPU_/(Nx_*Ny_*iter) << " (s)\n";
  cout<< left << setw(35) <<"\t Mega Lattice Site Update: "<< setw(8) << fixed      << (Nx_*Ny_*iter)/tCPU_/1000000.0 << " (MLSU/s)\n\n";

  // ======================= Simulation Parameters ======================= //
  
  writeSimulationParameters(iter);

  // ============================ Free memory ============================ //

  IntQuant_file.close();
  
  delete[] ex_; ex_   = 0;
  delete[] ey_; ey_   = 0;
  delete[] w_;  w_    = 0;
  delete[] Hxx_; Hxx_ = 0;
  delete[] Hxy_; Hxy_ = 0;
  delete[] Hyy_; Hyy_ = 0;

  deallocateDbl(Nx_, rho); deallocateDbl(Nx_, ux); deallocateDbl(Nx_, uy);
  deallocateDbl(Nx_, ux0); deallocateDbl(Nx_, uy0);
  deallocateDbl(Nx_, Ny_, f); deallocateDbl(Nx_, Ny_, fEq);
  deallocateDbl(Nx_, Ny_, fStream);
  // ++++++++++++++++++++++++++++++++++++++++++++++ //

  return 0;

}
