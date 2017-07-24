#ifndef INPUTFILE_HPP_INCLUDED
#define INPUTFILE_HPP_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace std;


void getMeshData(string& vtkFileName);

void initFromConfigFile(string& vtkFileName, double **rho, double **ux, double **ux0,
  					    double **uy, double **uy0);

void convertToLattice(double **rho, double **ux, double **uy, double **ux0, double **uy0);


#endif // INPUTFILE_HPP_INCLUDED
