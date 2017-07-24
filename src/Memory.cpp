#include "./Memory.hpp"
#include <iostream>
#include <cstring>

int ***allocateInt(int Nx, int Ny, int Q){
  int ***data  = new int**[Nx];
  for (int i = 0; i < Nx; ++i){
    data[i]  = new int*[Ny];
    for (int j = 0; j < Ny; ++j){
      data[i][j]  = new int[Q];
      std::memset(data[i][j], 0, sizeof(data[i][j]));  
    }
  }
  return data;
}

int **allocateInt(int Nx, int Ny){
  int **data  = new int*[Nx];
  for (int i = 0; i < Nx; ++i){
    data[i]  = new int[Ny];
    std::memset(data[i], 0, sizeof(data[i]));  
  }
  return data;
}

double ***allocateDbl(int Nx, int Ny, int Q){
  double ***data  = new double**[Nx];
  for (int i = 0; i < Nx; ++i){
    data[i]  = new double*[Ny];
    for (int j = 0; j < Ny; ++j){
      data[i][j]  = new double[Q];
      std::memset(data[i][j], 0, sizeof(data[i][j]));  
    }
  }
  return data;
}

double **allocateDbl(int Nx, int Ny){
  double **data  = new double*[Nx];
  for (int i = 0; i < Nx; ++i){
    data[i]  = new double[Ny];
    std::memset(data[i], 0, sizeof(data[i]));  
  }
  return data;
}

void deallocateInt(int Nx, int **t){
  for (int i=0; i < Nx; i++){
    delete[] t[i];
  }
  delete[] t;
  t = 0;
}

void deallocateInt(int Nx, int Ny, int ***t){
  for(int i(0); i < Nx; i++){
    for(int j(0); j < Ny; j++){
      delete[] t[i][j];
    }
    delete[] t[i];
  }
  delete[] t;
  t = 0;
}

void deallocateDbl(int Nx, double **t){
  for (int i=0; i < Nx; i++){
    delete[] t[i];
  }
  delete[] t;
  t = 0;
}

void deallocateDbl(int Nx, int Ny, double ***t){
  for(int i(0); i < Nx; i++){
    for(int j(0); j < Ny; j++){
      delete[] t[i][j];
    }
    delete[] t[i];
  }
  delete[] t;
  t = 0;
}

