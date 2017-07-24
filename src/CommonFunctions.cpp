#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "CommonFunctions.hpp"
#include "GlobalVars.hpp"
#include  <omp.h>

using namespace std;

// ========================== Lattice Generation =========================== //

void generateLattice(){

  ex_  = new double[Q_];
  ey_  = new double[Q_];
  w_   = new double[Q_];
  Hxx_ = new double[Q_];
  Hxy_ = new double[Q_];
  Hyy_ = new double[Q_];

  if (LBMmodel_ == "D2Q4"){
    ex_[0] = 0.; ey_[0] = 0.; w_[0] = 0.;                 //         2
    ex_[1] = 1.; ey_[1] = 0.; w_[1] = 1./4.;              //         *
    ex_[2] = 0.; ey_[2] = 1.; w_[2] = 1./4.;              //     3 * * * 1
    ex_[3] =-1.; ey_[3] = 0.; w_[3] = 1./4.;              //         *
    ex_[4] = 0.; ey_[4] =-1.; w_[4] = 1./4.;              //         4
    cs_ = 1./sqrt(2.);                                   
  }
  else if (LBMmodel_ == "D2Q5"){
    ex_[0] = 0.; ey_[0] = 0.; w_[0] = 1./3.;              //         2
    ex_[1] = 1.; ey_[1] = 0.; w_[1] = 1./6.;              //         *
    ex_[2] = 0.; ey_[2] = 1.; w_[2] = 1./6.;              //     3 * 0 * 1
    ex_[3] =-1.; ey_[3] = 0.; w_[3] = 1./6.;              //         *
    ex_[4] = 0.; ey_[4] =-1.; w_[4] = 1./6.;              //         4
    cs_ = 1./sqrt(3.);                                   
  }
  else if (LBMmodel_ == "D2Q9"){
    ex_[0] = 0.; ey_[0] = 0.; w_[0] = 4./9.;
    ex_[1] = 1.; ey_[1] = 0.; w_[1] = 1./9.;              //   4-----3-----2
    ex_[2] = 1.; ey_[2] = 1.; w_[2] = 1./36.;             //   | *   *   * |
    ex_[3] = 0.; ey_[3] = 1.; w_[3] = 1./9.;              //   |   * * *   |
    ex_[4] =-1.; ey_[4] = 1.; w_[4] = 1./36.;             //   5 * * 0 * * 1
    ex_[5] =-1.; ey_[5] = 0.; w_[5] = 1./9.;              //   |   * * *   |
    ex_[6] =-1.; ey_[6] =-1.; w_[6] = 1./36.;             //   | *   *   * |
    ex_[7] = 0.; ey_[7] =-1.; w_[7] = 1./9.;              //   6-----7-----8
    ex_[8] = 1.; ey_[8] =-1.; w_[8] = 1./36.;
    cs_ = 1./sqrt(3.);                                   

    // Definition of Hermite polynomials
    double cs2 = cs_ * cs_;
    for (int a(0); a < Q_; a++){
      Hxx_[a] = ex_[a] * ex_[a] - cs2;
      Hxy_[a] = ex_[a] * ey_[a];
      Hyy_[a] = ey_[a] * ey_[a] - cs2;
    }
  }
}

// ============================= Equilibrium DF ============================ //

void EquilibriumDF(double **rho, double **ux, double **uy, 
                   double ***fEq){

  clock_t tmp_tEq = clock();
  double cs2 = cs_ * cs_;
  double cs4 = cs2 * cs2;
  double cs6 = cs4 * cs2;

  if (aeroOrder_ == 0){
    // Diffusion equation
    for (int i(0); i < Nx_; i++){
      for (int j(0); j < Ny_; j++){
        for (int a(0); a < Q_; a++){
          fEq[i][j][a] = w_[a] * rho[i][j] ; 
        }
      }
    }
  }

  if (aeroOrder_ == 1){
    //Advection-diffusion with u=cst
    for (int i(0); i < Nx_; i++){
      for (int j(0); j < Ny_; j++){
        for (int a(0); a < Q_; a++){
          double uE =  (ex_[a] * ux[i][j] + ey_[a] * uy[i][j]) / cs2;
          fEq[i][j][a] = w_[a] * rho[i][j] * (1. + uE); 
        }
      }
    }
  }

  if (aeroOrder_ == 2){
    // Low compressible isothermal NS equations
#pragma omp parallel for
    for (int i=0; i < Nx_; i++){
      for (int j=0; j < Ny_; j++){
        double uSq = (ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j]) / cs2;
        for (int a=0; a < Q_; a++){
          double uE = (ex_[a] * ux[i][j] + ey_[a] * uy[i][j]) / cs2;
          fEq[i][j][a] = w_[a] * rho[i][j] * (1. + uE + uE*uE/2. - uSq/2.);
        }
      }
    }
  }  
  
  tmp_tEq = clock() - tmp_tEq;
  tEq_ += (double) tmp_tEq / CLOCKS_PER_SEC;
}

// =============== Collision step ==================================== //
void collisionBGK(double ***f, double ***fEq){


  clock_t tmp_tColl = clock();
  double omega = 1. /tau_;

  // Collision for aero distributions
#pragma omp parallel for
  for (int i=0; i< Nx_; i++){
    for (int j=0; j < Ny_; j++){
      for (int a=0; a < Q_; a++){
        f[i][j][a] = f[i][j][a] - omega * (f[i][j][a] - fEq[i][j][a]);
      }
    }
  }

  tmp_tColl = clock() - tmp_tColl;
  tColl_ += (double) tmp_tColl / CLOCKS_PER_SEC;
}

// ============================ Streaming step ============================= //
void streaming(double ***f, double ***fStream){


  clock_t tmp_tStream = clock();
  int in, ip, jn, jp;

  if (LBMmodel_=="D2Q4"){
#pragma omp parallel for
    for (int i=0; i< Nx_; i++){       
      for (int j=0; j < Ny_; j++){  
        ip = ((i+1+Nx_))%Nx_;
        in = ((i-1+Nx_))%Nx_;
        jp = ((j+1+Ny_))%Ny_;
        jn = ((j-1+Ny_))%Ny_;

        fStream[ip][j][1] = f[i][j][1];
        fStream[i][jp][2] = f[i][j][2];
        fStream[in][j][3] = f[i][j][3];
        fStream[i][jn][4] = f[i][j][4];

      }
    }
  }
      
  if (LBMmodel_=="D2Q5"){
#pragma omp parallel for
    for (int i=0; i< Nx_; i++){       
      for (int j=0; j < Ny_; j++){  
        ip = ((i+1+Nx_))%Nx_;
        in = ((i-1+Nx_))%Nx_;
        jp = ((j+1+Ny_))%Ny_;
        jn = ((j-1+Ny_))%Ny_;

        fStream[i][j][0] = f[i][j][0];
        fStream[ip][j][1] = f[i][j][1];
        fStream[i][jp][2] = f[i][j][2];
        fStream[in][j][3] = f[i][j][3];
        fStream[i][jn][4] = f[i][j][4];

      }
    }
  }

  if (LBMmodel_=="D2Q9"){
#pragma omp parallel for
    for (int i=0; i< Nx_; i++){       
      for (int j=0; j < Ny_; j++){  
        ip = ((i+1+Nx_))%Nx_;
        in = ((i-1+Nx_))%Nx_;
        jp = ((j+1+Ny_))%Ny_;
        jn = ((j-1+Ny_))%Ny_;

        fStream[i][j][0]   = f[i][j][0];
        fStream[ip][j][1]  = f[i][j][1];
        fStream[ip][jp][2] = f[i][j][2];
        fStream[i][jp][3]  = f[i][j][3];
        fStream[in][jp][4] = f[i][j][4];
        fStream[in][j][5]  = f[i][j][5];
        fStream[in][jn][6] = f[i][j][6];
        fStream[i][jn][7]  = f[i][j][7];
        fStream[ip][jn][8] = f[i][j][8];
      }
    }
  }

#pragma omp parallel for
  for (int i=0; i < Nx_; i++){
    for (int j=0; j < Ny_; j++){
      for (int a=0; a < Q_; a++){
        f[i][j][a] = fStream[i][j][a];
      }
    }
  }

  tmp_tStream = clock() - tmp_tStream;
  tStream_ += (double) tmp_tStream / CLOCKS_PER_SEC;
}

// ===================== Compute macroscopic variable ====================== //
void macroVar(double **rho, double **ux, double **uy, double ***f){


  clock_t tmp_tMacro = clock();

  if (aeroOrder_ == 0 || aeroOrder_ == 1){
    // Diffusion or Advection-Diffusion -> Computation of rho only !
#pragma omp parallel for
    for (int i=0; i < Nx_; i++){
      for (int j=0; j < Ny_; j++){
        rho[i][j] = 0.;
        for (int a=0; a < Q_; a++){
          rho[i][j] += f[i][j][a];
        }
      }
    }
  }

  if (aeroOrder_ == 2){
    // Isothermal low compressible NS -> Computation of rho and u
#pragma omp parallel for
    for (int i=0; i < Nx_; i++){
      for (int j=0; j < Ny_; j++){
        rho[i][j] = 0.;
        ux[i][j] = 0.;
        uy[i][j] = 0.;
        for (int a=0; a < Q_; a++){
          rho[i][j] += f[i][j][a];
          ux[i][j]  += ex_[a] * f[i][j][a];
          uy[i][j]  += ey_[a] * f[i][j][a];
        }
        ux[i][j] /= rho[i][j];       // we come back to velocity here !
        uy[i][j] /= rho[i][j];
      }
    }
  }

  tmp_tMacro = clock() - tmp_tMacro;
  tMacro_ += (double) tmp_tMacro / CLOCKS_PER_SEC;
}

// ======================== 2nd order regularization ======================= //
// Jonas Latt Regularization ("Lattice Boltzmann method with regularized
// pre-collision distribution functions", 2006)
void regularization(double ***f, double ***fEq){

  clock_t tmp_tReg = clock();

  // The aim here is to reconstruct f=fEq + fneq
  // 1) Computation of off-equilibrium coefficients a1xx, a1xy and a1yy
  //    by projection onto Hermite polynomials 
  // 2) Computation of fneq thanks to the Hermite polynomial development
  // 3) Reconstruction of f = fEq + fneq
  // Second order Hermite polynomials are already implemented as global variables in the code :
  // Hxx_[a], Hxy_[a], Hyy_[a] where a stands for the velocity index

  double cs2 = cs_ * cs_ ;
  double cs4 = cs2 * cs2 ;
  double fneq(0), a1xx(0), a1xy(0), a1yy(0);
  double omega = 1./tau_;

#pragma omp parallel for
  for (int i=0; i < Nx_; i++){
    for (int j=0; j < Ny_; j++){
     a1xx =0;
     a1xy =0;
     a1yy=0;
      for (int a=0; a < Q_; a++){
        a1xx += (f[i][j][a] -fEq[i][j][a])*Hxx_[a];
        a1xy += (f[i][j][a] -fEq[i][j][a])*Hxy_[a];
        a1yy += (f[i][j][a] -fEq[i][j][a])*Hyy_[a];
      }
      for (int a=0; a < Q_; a++){
        fneq = (w_[a]/(2*cs4))*(a1xx*Hxx_[a] +2*a1xy*Hxy_[a] +a1yy*Hyy_[a]);
        f[i][j][a] = fEq[i][j][a] +fneq;
      }
    }
  }

  tmp_tReg = clock() - tmp_tReg;
  tReg_ += (double) tmp_tReg / CLOCKS_PER_SEC;
}
