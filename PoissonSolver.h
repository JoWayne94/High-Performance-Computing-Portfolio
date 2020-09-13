/**
 * @file PoissonSolver.h
 *
 * High-performance Computing
 *
 * Header file for PoissonSolver.cpp
 *
 * Solves the Poisson problem with Cholesky factorization
 */
#ifndef POISSON_SOLVER
#define POISSON_SOLVER

class PoissonSolver {
public:
	PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy);
	~PoissonSolver(); 
	
	void BoundariesPartition(const double* top, const double* bottom, const double* left, const double* right);
	void SolvePoisson(double* x, const double* f);
	
private:
	double* A = nullptr;	
	double* b = nullptr;  // boundary vector for cblas routines
	int ld_A;
	int Nx;
	int Ny;
	double dx;
	double dy;
	int n_dim;
	int info;
};

#endif































