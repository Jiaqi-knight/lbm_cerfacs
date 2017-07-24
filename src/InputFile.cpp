#include <iostream>
#include <cmath>
// #include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <sys/stat.h>
#include "InputFile.hpp"
#include "GlobalVars.hpp"

using namespace std;

void getMeshData(string& vtkFileName){

  ifstream inputFile;
  inputFile.open(&vtkFileName[0]);

  if (inputFile.fail()){
    cout << "Error : vtk file doesn't exist !" << endl;
    exit (EXIT_FAILURE);
  }

  string text;
  inputFile >> text;

  while(text!="DIMENSIONS"){
    getline(inputFile, text);
    inputFile >> text;
  }
  inputFile >> text;
  Nx_ = atoi(text.c_str());
  inputFile >> text;
  Ny_ = atoi(text.c_str());

  while(text!="SPACING"){
    getline(inputFile, text);
    inputFile >> text;
  }
  inputFile >> text;
  dx_ = atof(text.c_str());

  inputFile.close();
}


void initFromConfigFile(string& vtkFileName, double **rho, double **ux, double **ux0,
                        double **uy, double **uy0){

  ifstream inputFile;
  inputFile.open(&vtkFileName[0]);
  double tmp_val;

  string text;
  inputFile >> text;

  // Now read in parameters
  while(inputFile){
    while(text!="SCALARS"){
      getline(inputFile, text);
      inputFile >> text;
    }
    inputFile >> text;

    // Read pressure
    double pressure;
    if (text == "Pressure"){
      getline(inputFile, text);
      getline(inputFile, text);
      for (int j(0); j < Ny_; j++){
        for (int i(0); i < Nx_; i++){
          inputFile >> pressure;
          rho[i][j] = pressure / (Rgas_ * Tref_) ; 
          // Density/Pressure and Temperature are decoupled in the present version !
        }
      }
      inputFile >> text;
    }

    // Read (dimensional) density
    if (text == "Density"){
      getline(inputFile, text);
      getline(inputFile, text);
      for (int j(0); j < Ny_; j++){
        for (int i(0); i < Nx_; i++){
          inputFile >> rho[i][j];
        }
      }
      inputFile >> text;
    }

    // Read velocity_X
    if (text == "velocity_X"){
      getline(inputFile, text);
      getline(inputFile, text);
      for (int j(0); j < Ny_; j++){
        for (int i(0); i < Nx_; i++){
          inputFile >> ux[i][j];
          ux0[i][j] = ux[i][j];
        }
      }
      inputFile >> text;
    }
  
    // Read velocity_Y
    if (text == "velocity_Y"){
      getline(inputFile, text);
      getline(inputFile, text);
      for (int j(0); j < Ny_; j++){
        for (int i(0); i < Nx_; i++){
          inputFile >> uy[i][j];
          uy0[i][j] = uy[i][j];
        }
      }
      inputFile >> text;
    }
  }
  inputFile.close();
}

void convertToLattice(double **rho, double **ux, double **uy, double **ux0, double **uy0){
  
  //GenerateLattice has to be computed first to have access to correct value of cs_ 
  lx_ = (Nx_-1) * dx_;
  ly_ = (Ny_-1) * dx_;

  dt_ = cryoRatio_ * cs_ * dx_ / sqrt(Gamma_ * Rgas_ * Tref_);
  // time step computed thanks to the following formula : 
  // C = cryoRatio * cs_ * dx_ / dt_
  // C : physical sound speed, cs_ : lattice sound speed (1./sqrt(3.) for D2Q9)

  if (simTime_ > 0.) nite_ = simTime_ / dt_ ;

  NuLb_ = Nu_ * dt_ / (dx_ * dx_);             // Kinematic viscosity in lattice units
  // NuLb_ is the non-dimensional viscosity

  tau_ = NuLb_/(cs_ * cs_) + 0.5;                //  Non-dimensional relaxation time for D2Q9
  // Demonstration : mu = rho * cs_ * cs_ * tau_ (definition of viscosity)
  // => tau_ = NuLb_ / (cs_ * cs_)
  // The term "+ 0.5" appears after the explicit variable change


  // Now convert macroscopic quantities into lattice units
  for (int i(0); i < Nx_; i++){
    for (int j(0); j < Ny_; j++){
      rho[i][j] *= Rgas_ * Tref_ / Pref_;   // OK for decoupled lattices.
      ux[i][j]  *= dt_ / dx_;
      uy[i][j]  *= dt_ / dx_;
      ux0[i][j]  = ux[i][j];
      uy0[i][j]  = uy[i][j];
    }
  }
}
