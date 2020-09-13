/**
 * @file LidDrivenCavity.h
 *
 * High-Performance Computing
 *
 * Header file of LidDrivenCavity.cpp
 *
 * Defines all public classes, private variables and the poisson solver object
 *
 */
#ifndef LID_DRIVEN_CAVITY
#define LID_DRIVEN_CAVITY

#include "PoissonSolver.h"

#include <string>
#include <mpi.h>

using namespace std;

class LidDrivenCavity {
public:
    LidDrivenCavity(MPI_Comm gridcomm, int rank, int* cart_coords, double* global_coords, int* neighbour_rank, bool &dt_bool, double deltat, double finalt, int nx, int ny, double re, double dx, double dy);
    ~LidDrivenCavity();
    
    void SetDomainSize(double xlen, double ylen);
    void Initialise();
    void Integrate();

    // Solver
    void CreateBandedMatrices();
    void BoundariesPartitionA(double* b, char x);
    void BoundariesPartitionB(double* b, char x);
    void BoundariesPartitionC(double* b, char x);
    void CalculateVelocity();
    void VorticityBCs();
    void InteriorVorticity();
    void NextInteriorVorticity();
    void Algorithm();

    // MPI
    void SendRecvProcess(double* x, double* x_top, double* x_bottom, double* x_left, double* x_right);
    void Output(int Px, int Py, double Lx, double Ly);

private:
    double* v = nullptr;
    double* s = nullptr;

    double dt;
    double T;
    unsigned int Nx;
    unsigned int Ny;
    double Lx;
    double Ly;
    double Re;

    double* A = nullptr;	
    double* B = nullptr;    
    double* C = nullptr;
    int n_dim;
    int ld_A; 
    int ld_B; 
    int ld_C;

    double* vbc_top = nullptr;
    double* vbc_bottom = nullptr;
    double* vbc_left = nullptr;
    double* vbc_right = nullptr;
    double* sbc_top = nullptr;
    double* sbc_bottom = nullptr;
    double* sbc_left = nullptr;
    double* sbc_right = nullptr;

    double* s_prev = nullptr;
    double* Ux = nullptr;
    double* Uy = nullptr;

    double dx;
    double dy;
    const double U = 1.0;  			// Lid velocity is 1
    PoissonSolver *Poisson_Solver = nullptr;	// poisson solver object

    // MPI configuration 
    MPI_Comm gridcomm;
    int rank;
    int cart_coords[2];       // Coordinates in cartesian topology
    double global_coords[2];  // Global coordinates
    int neighbour_rank[4];    // neighbour_rank ranks
};

#endif