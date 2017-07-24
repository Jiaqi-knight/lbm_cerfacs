#ifndef MEMORY_HPP_INCLUDED
#define MEMORY_HPP_INCLUDED


int ***allocateInt(int Nx, int Ny, int Q);

int **allocateInt(int Nx, int Ny);

double ***allocateDbl(int Nx, int Ny, int Q);

double **allocateDbl(int Nx, int Ny);

void deallocateInt(int Nx, int **t);

void deallocateInt(int Nx, int Ny, int ***t);

void deallocateDbl(int Nx, double **t);

void deallocateDbl(int Nx, int Ny, double ***t);


#endif // MEMORY_HPP_INCLUDED
