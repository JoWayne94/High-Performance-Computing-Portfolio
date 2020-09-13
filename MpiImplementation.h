/**
 * @file MpiImplementation.h
 *
 * High-Performance Computing
 *
 * Header file of MpiImplementation.cpp
 *
 * Declare all functions used to control MPI functionalities
 *
 */
#ifndef MPI_CONFIG
#define MPI_CONFIG

#include <mpi.h>

bool Processors(const int &np, const int &Px, const int &Py);
void Neighbours(MPI_Comm gridcomm, int* neighbour_rank);
void Allocation(const int &rank, const int &Nx, const int &Ny, const int &Px, const int &Py, const double &dx, const double &dy, int* cart_coords, double* global_coords, int &nx, int &ny);

#endif