#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include "Memory.hpp"
#include "Output.hpp"
#include "GlobalVars.hpp"

using namespace std;


void writeOutputVTK(int iter, double **rho, double **ux, double **uy){

  clock_t tmp_tOut = clock();
  double uNorm, Mach;
  double **vortZ = allocateDbl(Nx_, Ny_);
  computeVorticity(ux, uy, vortZ);

  string text("\n  Writing initial results...\n");

  if (iter != 0) text = " || Writing results...";
  cout << text;
  mkdir("output", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::stringstream fNameStream;
  double eps = 0.0000000000000001;
  fNameStream << "output/sol_" << std::setfill('0') << std::setw(7) << iter << ".vtk";
  string fName = fNameStream.str();
  ofstream resultFlux;
  resultFlux.open(fName.c_str());
  resultFlux << "# vtk DataFile Version 2.0" << endl ;
  resultFlux << "iteration " << iter << endl ;
  resultFlux << "ASCII" << endl ;
  resultFlux << endl ;
  resultFlux << "DATASET STRUCTURED_POINTS" << endl ;
  resultFlux << "DIMENSIONS " << Nx_ << " " << Ny_ << " 1" << endl ;
  resultFlux << "ORIGIN 0 0 0" << endl ;
  resultFlux << "SPACING " << dx_ << " " << dx_ << " " << dx_ << endl ;
  resultFlux << endl ;
  resultFlux << "POINT_DATA " << Nx_*Ny_ << endl ;

  // Density
  resultFlux << "SCALARS Density float 1" << endl ;
  resultFlux << "LOOKUP_TABLE default" << endl ;
  for (int j(0); j < Ny_; j++){
    for (int i(0); i < Nx_; i++){
      if (abs(rho[i][j]) > eps){ // Fix for Paraview issue when value very low
        resultFlux << rho[i][j] * Pref_ / (Rgas_*Tref_) << " ";
      }
      else{
        resultFlux << 0.0 << " ";
      }
    }
      resultFlux << endl ;
  }
  resultFlux << endl ;

  // Pressure
  resultFlux << "SCALARS Pressure float 1" << endl ;
  resultFlux << "LOOKUP_TABLE default" << endl ;
  for (int j(0); j < Ny_; j++){
    for (int i(0); i < Nx_; i++){
      if (abs(rho[i][j]) > eps){ // Fix for Paraview issue when value very low
        // A modifier
        resultFlux << rho[i][j]*Pref_ << " ";
      }
      else{
        resultFlux << 0.0 << " ";
      }
    }
      resultFlux << endl ;
  }
  resultFlux << endl ;

  // velocity_X
  resultFlux << "SCALARS velocity_X float 1" << endl ;
  resultFlux << "LOOKUP_TABLE default" << endl ;
  for (int j(0); j<Ny_; j++){
    for (int i(0); i<Nx_; i++){
      if (abs(ux[i][j]) > eps){ // Fix for Paraview issue when value very low
        resultFlux << ux[i][j] * dx_ / dt_ << " ";
      }
      else{
          resultFlux << 0.0 << " ";
      }
    }
    resultFlux << endl ;
  }
  resultFlux << endl ;

  // velocity_Y
  resultFlux << "SCALARS velocity_Y float 1" << endl ;
  resultFlux << "LOOKUP_TABLE default" << endl ;
  for (int j(0); j<Ny_; j++){
    for (int i(0); i<Nx_; i++){
      if (abs(uy[i][j]) > eps){ // Fix for Paraview issue when value very low
        resultFlux << uy[i][j] * dx_ / dt_ << " ";
      }
      else{
          resultFlux << 0.0 << " ";
      }
    }
    resultFlux << endl ;
  }
  resultFlux << endl ;

  // Mach number
  resultFlux << "SCALARS Mach float 1" << endl ;
  resultFlux << "LOOKUP_TABLE default" << endl ;
  for (int j(0); j<Ny_; j++){
    for (int i(0); i<Nx_; i++){
      uNorm = sqrt(ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j]);
      Mach = uNorm / cs_;
      if (abs(Mach) > eps){ // Fix for Paraview issue when value very low
        resultFlux << Mach << " ";
      }
      else{
          resultFlux << 0.0 << " ";
      }
    }
    resultFlux << endl ;
  }
  resultFlux << endl ;

  // vorticity_Z
  resultFlux << "SCALARS vorticity_Z float 1" << endl ;
  resultFlux << "LOOKUP_TABLE default" << endl ;
  for (int j(0); j<Ny_; j++){
    for (int i(0); i<Nx_; i++){

      if (abs(vortZ[i][j]) > eps){ // Fix for Paraview issue when value very low
        resultFlux << vortZ[i][j] / dt_ << " ";
      }
      else{
          resultFlux << 0.0 << " ";
      }
    }
    resultFlux << endl ;
  }
  resultFlux << endl ;

  tmp_tOut = clock() - tmp_tOut;
  tOut_ += (double) tmp_tOut / CLOCKS_PER_SEC;
}

void writeSimulationParameters(int t){
  
  string const dataFile("./output/SimulationParameters.log");
  ofstream myFlux(dataFile.c_str()); // we try to open the data file 
                                     // (c_str() is used to get the name)
  if (LBMmodel_=="D2Q4"){
  myFlux << "//=========== D2Q4 Computation =============//" << endl << endl ;
  }
  if (LBMmodel_=="D2Q5"){
  myFlux << "//=========== D2Q5 Computation =============//" << endl << endl ;
  }
  if (LBMmodel_=="D2Q9"){
  myFlux << "//========== D2Q9 Computation ============//" << endl << endl ;
  }
  
  double c  = cs_ * dx_/dt_;
  myFlux << "//========== Physical Parameters ===========//" << endl << endl;
  myFlux << "X Length :                    " << lx_ << " m" << endl;
  myFlux << "Y Length :                    " << ly_ << " m" << endl;
  myFlux << "Reference Temperature :       " << Tref_ << " K" << endl;
  myFlux << "Reference Pressure :          " << Pref_ << " Pa" << endl;
  myFlux << "Sound Speed :                 " << c << " m/s" << endl;
 // myFlux << "Mean Velocity (m/s):          " << ly << endl;

  myFlux << "Kinematic Viscosity :         " << Nu_ << " m2/s" << endl;
  //myFlux << "Reynolds Number :             " << *Re << endl;
  //myFlux << "Mach Number :                 " << Ma << endl;
  //myFlux << "Reference Velocity :          " << *Uo << endl;

  myFlux << "//========== Numerical Parameters ===========//" << endl << endl;
  myFlux << "Space Step:                   " << dx_ << " m" << endl;
  myFlux << "Time  Step:                   " << dt_ << " s" << endl;
  myFlux << "Relaxation Time (Tau) :       " << tau_ << endl << endl << endl;

  double N = Nx_ * Ny_;

  myFlux << "//============ CPU Performances =============//" << endl << endl;
  myFlux << "Number of nodes :             " << N << endl;
  myFlux << "Number of iterations :        " << t << " ite" << endl;
  myFlux << "Mean time:           " << endl;
  myFlux << "    - Total :                 " << tCPU_    << " (s)" << endl;
  myFlux << "    - Equilibrium :           " << tEq_     << " (s)" << endl;
  myFlux << "    - Regularization :        " << tReg_    << " (s)" << endl;
  myFlux << "    - Collision :             " << tColl_   << " (s)" << endl;
  myFlux << "    - Streaming :             " << tStream_ << " (s)" << endl;
  myFlux << "    - Macroscopic :           " << tMacro_  << " (s)" << endl;
  myFlux << "    - Output :                " << tOut_    << " (s)" << endl;
}
void printHeader(){
  cout << "       *************************************************************" << endl;
  cout << "       *                                                           *" << endl;
  cout << "       *     __      ______  ____    ____    ___     ____    ____  *" << endl;
  cout << "       *    / /     / ____/ / __ \\  / __ \\  /   |   / __ \\  / __ \\ *" << endl;
  cout << "       *   / /     / __/   / / / / / /_/ / / /| |  / /_/ / / / / / *" << endl;
  cout << "       *  / /___  / /___  / /_/ / / ____/ / ___ | / _, _/ / /_/ /  *" << endl;
  cout << "       * /_____/ /_____/  \\____/ /_/     /_/  |_|/_/ |_| /_____/   *" << endl;
  cout << "       *                                                           *" << endl;
  cout << "       *         \033[0m \033[1mL\033[0mattic\033[0m\033[1mE\033[0m b\033[0m\033[1mO\033[0mltzmann \033[0m \033[1mP\033[0ml\033[0m\033[1mA\033[0mtefo\033[0m\033[1mR\033[0mm \033[0m\033[1mD\033[0mevelopment         *" << endl;
  cout << "       *************************************************************" << endl;
  cout << "       *                   Copyright 2015, CERFACS                 *" << endl;
  cout << "       *************************************************************\033[0m" << endl;
}

void printInfo(){
  cout << "\n";
  cout << "       ============================================================="<< endl;
  cout << "       =                        Read parameters                    ="<< endl;
  cout << "       ============================================================="<< endl;
  cout << left << setw(32) << "\t\tLBM modelisation: "             << setw(7) << LBMmodel_          << endl;
  cout << left << setw(32) << "\t\tRegularisation: "               << setw(7) << BoolToString(reg_) << endl;
  cout << left << setw(32) << "\t\tReference Temperature: "        << setw(7) << Tref_              << endl;
  cout << left << setw(32) << "\t\tReference Pressure: "           << setw(7) << Pref_              << endl;
  cout << left << setw(32) << "\t\tKinematic Viscosity: "          << setw(7) << Nu_                << endl;
  if (nite_ > 0){
    cout << left << setw(32) << "\t\tNumber of iterations:"          << setw(7) << nite_              << endl;
  } else if (simTime_ > 0.){
    cout << left << setw(32) << "\t\tSimulation time (s):"           << setw(7) << simTime_           << endl;
  }
  
  cout << left << setw(32) << "\t\tResult output frequency: "      << setw(7) << freqOut_           << endl;
  cout << left << setw(32) << "\t\tConfiguration file: "           << setw(7) << config_file        << endl;
  cout << left << setw(32) << "\t\tR gas constant: "               << setw(7) << Rgas_              << endl;
  cout << left << setw(32) << "\t\tGamma: "                        << setw(7) << Gamma_             << endl;
  cout << left << setw(32) << "\t\tAero distribution dev. order: " << setw(7) << aeroOrder_         << endl;
  cout << "       ============================================================="<< endl;
}

void computeVorticity(double **ux, double **uy, double **vortZ){
  for (int j(0); j<Ny_; j++){
    for (int i(0); i<Nx_; i++){
      int in = (i-1+Nx_)%Nx_;
      int ip = (i+1+Nx_)%Nx_;
      int jn = (j-1+Ny_)%Ny_;
      int jp = (j+1+Ny_)%Ny_;

      double dudy = (ux[i][jp] - ux[i][jn]) / 2.;
      double dvdx = (uy[ip][j] - uy[in][j]) / 2.;

      vortZ[i][j] = dvdx - dudy; 
    }
  }
}

void computeIntegratedQuantities(ofstream& IntQuant_file, double **rho, double **ux, double **uy, double time){

  if (time == 0.){
    IntQuant_file << "#time meanKineticEnergy RMSKineticEnergy meanEnstrophy RMSEnstrophy " << endl;
  }

  double **vortZ = allocateDbl(Nx_, Ny_);
  computeVorticity(ux, uy, vortZ);

  double meanKineticEnergy = 0.;
  double meanEnstrophy = 0.;
  for (int j(0); j<Ny_; j++){
    for (int i(0); i<Nx_; i++){
      meanKineticEnergy += 0.5*rho[i][j]*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]);
      meanEnstrophy += 0.5 * vortZ[i][j] * vortZ[i][j];
    }
  }
  meanKineticEnergy /= (Nx_*Ny_);
  meanEnstrophy /= (Nx_*Ny_);

  double RMSKineticEnergy = 0.;
  double RMSEnstrophy = 0.;
  for (int i(0); i < Nx_; i++){
    for (int j(0); j < Ny_; j++){
      RMSKineticEnergy += (0.5*rho[i][j]*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]) - meanKineticEnergy)
                          * (0.5*rho[i][j]*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]) - meanKineticEnergy);
      RMSEnstrophy += (0.5*vortZ[i][j]*vortZ[i][j] - meanEnstrophy) 
                      * (0.5*vortZ[i][j]*vortZ[i][j] - meanEnstrophy);
    }
  }
  RMSKineticEnergy /= (Nx_*Ny_);
  RMSKineticEnergy = sqrt(RMSKineticEnergy);
  RMSEnstrophy /= (Nx_*Ny_); 
  RMSEnstrophy = sqrt(RMSEnstrophy);

  // We come back to dimensional quantities here
  meanKineticEnergy *= dx_ * dx_ / (dt_ * dt_) * Pref_ / (Rgas_*Tref_);
  RMSKineticEnergy  *= dx_ * dx_ / (dt_ * dt_) * Pref_ / (Rgas_*Tref_);
  meanEnstrophy     *= dx_ * dx_ / (dt_ * dt_);
  RMSEnstrophy      *= dx_ * dx_ / (dt_ * dt_);


  IntQuant_file << time << " " << meanKineticEnergy << " " << RMSKineticEnergy 
                  << " " << meanEnstrophy << " " << RMSEnstrophy <<  endl;

}
