/**
 * @file PoissonSolver.cpp
 *
 * High-performance Computing
 *
 * Solution of a Poisson problem to compute the stream-function at time t + ∆t
 *
 * Solves system of equations using a linear solver, cblas and lapack routines
 */
#include "PoissonSolver.h"
#include "cblas.h"

#include<cstring>
#include<iostream>
#include<iomanip>
#include<cmath>

using namespace std;

// Calls LaPack
#define F77NAME(x) x##_
extern "C" {
	// Computes the Cholesky factorization of a real symmetric positive definite band matrix A
	void F77NAME(dpbtrf) (const char &UPLO, const int &N, const int &KD,
			     const double *AB, const int &LDAB,	int &INFO);

	// Solves a system of linear equations A*X = B with a symmetric 
	// positive definite band matrix A using the Cholesky factorization
	void F77NAME(dpbtrs) (const char &UPLO, const int &N, const int &KD,
			     const int &NRHS, const double *AB, const int &LDAB,
			     double *B, const int &LDB, int &INFO);
}

/**
 * @brief Constructor
 *
 * @param   Nx      Number of grid points in x-direction
 * @param   Ny      Number of grid points in y-direction
 * @param   dx      Grid spacing in x-direction
 * @param   dy      Grid spacing in y-direction
 */
PoissonSolver::PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy) {
	this-> Nx = Nx;
	this-> Ny = Ny;
	this-> dx = dx;
	this-> dy = dy;

	// Create banded matrix A of the system
	n_dim = Nx * Ny;
	ld_A = Ny + 1;
	A = new double[n_dim * ld_A];
	cblas_dscal(n_dim * ld_A, 0.0, A, 1);

	for (unsigned int j = 0; j < n_dim; j++) {
		A[ld_A - 1 + j * ld_A] 		    	  = 2/(dx*dx) + 2/(dy*dy);
		if (j % Ny != 0)   A[ld_A - 2 + j * ld_A] = -1/(dy*dy);
		if (j >= Ny)       A[j * ld_A] 	     	  = -1/(dx*dx);
	}
	info = 1;
}

/**
 * @brief Destructor
 *
 * 	  Clean up memory
 */
PoissonSolver::~PoissonSolver() {
	delete[] A;
	delete[] b;
}

/**
 * @brief Accept boundary conditions from 4 walls and constructs boundary vector b
 *
 * @param   top    	Top boundary array of a partition
 * @param   bottom    	Bottom boundary array of a partition
 * @param   left   	Left most boundary array of a partition
 * @param   right  	Right most boundary array of a partition
 */
void PoissonSolver::BoundariesPartition(const double* top, const double* bottom, const double* left, const double* right) {
	b = new double[n_dim];
	cblas_dscal(n_dim, 0.0, b, 1);

	for (unsigned int i = 0; i < Nx; i++) {
		b[i * Ny] += bottom[i] * -1/(dy*dy);		// b vector inputs for bottom wall 
		b[i * Ny + Ny - 1] += top[i] * -1/(dy*dy);      // b vector inputs for top wall
	}

	for (unsigned int i = 0; i < Ny; i++) {
		b[i] += left[i] * -1/(dx*dx); 		        // b vector inputs for left wall
		b[i + (Nx - 1) * Ny] += right[i] * -1/(dx*dx);  // b vector inputs for right wall
	}
}

/**
 * @brief Solve the matrix problem \f$ Ax = b \f$ using the LAPACK routine DPBTRF & DPBTRS
 * 
 * @param   x		Update stream function, hence parameter x is s in main code
 * @param   bot    	Bottom boundary array of a partition
 */
void PoissonSolver::SolvePoisson(double* x, const double* f) {

	// Ensure factorisation occurs only once
	if (info != 0) {
		F77NAME(dpbtrf) ('U', n_dim , Ny, A, ld_A, info);
	}

	// Set RHS vector
	cblas_daxpy(n_dim, -1.0, f, 1, b, 1);
	cblas_dscal(n_dim, -1.0, b, 1);

	// Solve Ax = RHS, output is b
	F77NAME(dpbtrs) ('U', n_dim, Ny, 1, A, ld_A, b, n_dim, info);

	// Outputs x
	cblas_dcopy(n_dim, b, 1, x, 1);
}
