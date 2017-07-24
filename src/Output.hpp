#ifndef OUTPUT_HPP_INCLUDED
#define OUTPUT_HPP_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

void writeOutputVTK(int iter, double **rho, double **ux, double **uy);

void writeSimulationParameters(int iter);

void printHeader();

void printInfo();

void computeVorticity(double **ux, double **uy, double **vortZ);

void computeIntegratedQuantities(ofstream& IntQuant_file, double **rho, double **ux, double **uy, double time);

inline const char * const BoolToString(bool b)
{
  return b ? "True" : "False";
};

#endif // OUTPUT_HPP_INCLUDED
