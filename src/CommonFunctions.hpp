#ifndef COMMON_HPP_INCLUDED
#define COMMON_HPP_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>


void generateLattice();

void EquilibriumDF(double **rho, double **ux, double **uy, double ***fEq);

void macroVar(double **rho, double **ux, double **uy, double ***f);

void collisionBGK(double ***f, double ***fEq);

void streaming(double ***f, double ***fStream);

void regularization(double ***f, double ***fEq);

#endif // COMMON_HPP_INCLUDED
