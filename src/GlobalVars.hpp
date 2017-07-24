#ifndef GLOBALVARS_HPP_INCLUDED
#define GLOBALVARS_HPP_INCLUDED

  using namespace std;

  extern string LBMmodel_;        // choice of code : D2Q4, D2Q5 or D2Q9
  extern int D_;
  extern int Q_;                  // lattice number of speeds
  extern double *ex_, *ey_;       // lattice speed coordinates
  extern double *w_;              // lattice weights
  extern double cs_;              // lattice sound speed;

  extern double lx_, ly_;         // Size of the computational domain in meters
  extern double dx_;              // Mesh size
  extern double dt_;              // time step
  extern int Nx_, Ny_;            // Number of cells
  extern int nite_;               // Number of time steps (nite=0: conv criteria)
  extern double simTime_;         // Simulation time (in seconds)

  extern double Nu_;              // Kinematic viscosity in physical units
  extern double NuLb_;            // Kinematic viscosity in lattice units
  extern double tau_;             // Relaxation time (classical LBM)

  extern double Tref_;            // Reference temperature
  extern double Pref_;            // Reference pressure

  extern int freqOut_;            // Output frequency
  extern int freqIntQuant_;       // frequency of integrated quantities

  extern int aeroOrder_;          // Development order of fEq
  extern double Rgas_;            // Ideal gas constant 
  extern double Gamma_;           // heat capacity ratio
  extern double cryoRatio_;       // Cryogenic ratio = Mach cryo / real Mach

  extern int nThreads_;     //number of threads openMP

  extern bool reg_;
  extern string config_file;

  extern double *Hxx_, *Hxy_, *Hyy_;

  extern double tCPU_; 
  extern double tEq_; 
  extern double tReg_; 
  extern double tColl_; 
  extern double tStream_; 
  extern double tMacro_; 
  extern double tOut_;

#endif // GLOBALVARS_HPP_INCLUDED