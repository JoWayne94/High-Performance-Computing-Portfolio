/**
 * @file PoissonSolver.h
 *
 * High-performance Computing
 *
 * Header file for PoissonSolver.cpp
 *
 * Defines the Poisson solver object
 */
#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H
#pragma once
#include <string>
using namespace std;

/*  For the 2D solver, we assume a problem of [A]{x} = {b}. Since the boundary conditions
 *  are known, we reduce this to solving for the interior points only.
 *  A second order central difference scheme is used for the spatial derivatives.
 *  The reduced problem becomes [AInt]{xInt} = {bInt}, applying known boundary conditions on {bInt}
 */

class PoissonSolver {

public:
	PoissonSolver();
	~PoissonSolver(); 
	
    	void GenerateScaLaPackAIntMatrix();
    	void GenerateLaPackAIntMatrix();
    	double* GetScaLaPackAIntMatrix();
    	int GetScaLaPackAIntMatrixNx();
    	int GetScaLaPackAIntMatrixNy();
    
    	void SetScaLaPackAIntMatrix(double* aint, int aintnx, int aintny);
    	void SetVariables(int nx, int ny, double diagVar, double firstdiagVar, double seconddiagVar);
    	void SetVectors(double* xint, double* bVec);
    	void Updatex(double* xVec);
    
    	// Pre-factor and solve functions
    	void InitialiseScaLaPack(int px, int py);
    	void PrefactorAIntMatrixParallel();
    	void SolveParallel();
    	void PrefactorAIntMatrixSerial();
    	void SolveSerial();
	
private:
    	// [A]{x} = {b} matrices and vectors
    	double* AInt = nullptr;
    	double* xInt = nullptr;
    	double* bInt = nullptr;
    
    	// Main and super diagonals of A
   	double diag;
    	double firstdiag;
    	double seconddiag;
    
    	int Nx;
    	int Ny;
    	int bIntNx;
    	int bIntNy;
    	int AIntNx;
    	int AIntNy;
    
    	// MPI variables
    	int rank;
    	int nprocs;
    
    	int Px;
    	int Py;
    	int mype;
    	int npe;
    	int ctx;
    	int nrow;
    	int ncol;
    	int myrow;
    	int mycol;
    
    	int m;
    	int n;
    	int nb;
    	int bwl;
    	int bwu;
    	int ja;
    	int desca[7];
    	int lda;

    	double* prefactoredAInt = nullptr;
    	int* ipiv = nullptr;
    	double* af = nullptr;
    	int laf;
};

#endif
