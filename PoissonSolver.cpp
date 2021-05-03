/**
 * @file PoissonSolver.cpp
 *
 * High-performance Computing
 *
 * Solution to the Poisson problem to compute the stream-function at time t + ∆t
 * Define class member functions for the PoissonSolver class and implements a serial and parallel solver for a 2D Poisson Equation (Dirichlet Problem)
 * Solves system of equations using a linear solver, Cblacs, LaPack, and ScaLaPack routines
 */
#include "PoissonSolver.h"
#include "cblas.h"
#include<mpi.h>
#include<cstring>
#include<iostream>
#include<iomanip>
#include<cmath>

using namespace std;

// Call fortran routines to be used in C++
#define F77NAME(x) x##_
extern "C" {
    // Performs LU factorisation of a general banded matrix in serial
    void F77NAME(dgbtrf) (const int& m, const int& n, const int& kl,
                          const int& ku, double* A, const int& lda,
                          int* ipiv, int& info);
    
    // Solves prefactored-system of equations in serial
    void F77NAME(dgbtrs) (const char& trans, const int& n, const int& kl,
                          const int &ku, const int& nrhs, const double* A,
                          const int& lda, const int* ipiv, double* b,
                          const int& ldb, int& info);

    // Performs LU factorisation of a general banded matrix in parallel
    void F77NAME(pdgbtrf)(const int& n, const int& bwl, const int& bwu,
                          double* A, const int& ja, int* desca,
                          int* ipiv, double* af, int& laf,
                          double* work, int& lwork, int& info);
    
    // Solves prefactored-system of equations in parallel
    void F77NAME(pdgbtrs)(const char& trans, const int& n, const int& bwl,
                          const int& bwu, const int& nrhs, double* A,
                          const int& ja, int* desca, int* ipiv,
                          double* b, const int& ib, int* descb, double* af, int& laf,
                          double* work, int& lwork, int& info);
    
    // Sends a matrix from current process to destination process
    void F77NAME(dgesd2d)(const int& icontxt, const int& m, const int& n,
                          double* A, const int& lda, const int& rdest,
                          int& cdest);
    
    // Receives a matrix from the source process to the current process
    void F77NAME(dgerv2d)(const int& icontxt, const int& m, const int& n,
                          double* A, const int& lda, const int& rsrc,
                          int& csrc);

    // Return number of available processes ready to use
    void Cblacs_pinfo(int*, int*);

    // Get values that BLACS uses for internal defaults
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
    void Cblacs_gridexit(int);
}

/**
 * @brief Constructor of the PoissonSolver class
 */
PoissonSolver::PoissonSolver() {
}

/**
 * @brief Destructor of the PoissonSolver class
 *
 * 	  Clean up memory
 */
PoissonSolver::~PoissonSolver() {

    // Deallocate memory of class arrays
    delete[] AInt;
    delete[] xInt;
    delete[] bInt;
    delete[] prefactoredAInt;
    delete[] ipiv;
    delete[] af;
    
    // Finalise cblacs
    Cblacs_gridexit(ctx);
}

/**
 * @brief Set variables
 *        The Int variables are reduced variables in the interior matrix
 *
 * @param   nx    		Grid spacing in the x-direction
 * @param   ny    		Grid spacing in the y-direction
 * @param   diagVar   		Main diagonal of A
 * @param   firstdiagVar  	Super diagonal of A, 1/dy/dy
 * @param   seconddiagVar  	Super diagonal of A, 1/dx/dx
 */
void PoissonSolver::SetVariables(int nx, int ny, double diagVar, double firstdiagVar, double seconddiagVar) {

    Nx = nx;
    Ny = ny;
    diag = diagVar;
    firstdiag = firstdiagVar;
    seconddiag = seconddiagVar;
    bIntNx = Nx - 2;
    bIntNy = Ny - 2;
    
    // Initialize interior arrays
    bInt = new double[bIntNx * bIntNy]{};
    xInt = new double[bIntNx * bIntNy]{};
}

/**
 * @brief Get generated AInt matrix
 *
 * @return  Pointer to AInt matrix
 */
double* PoissonSolver::GetScaLaPackAIntMatrix() {
    return AInt;
}

/**
 * @brief Get number of grid points in x-direction
 *
 * @return  Grid points in x-direction
 */
int PoissonSolver::GetScaLaPackAIntMatrixNx() {
    return AIntNx;
}

/**
 * @brief Get number of grid points in y-direction
 *
 * @return  Grid points in y-direction
 */
int PoissonSolver::GetScaLaPackAIntMatrixNy() {
    return AIntNy;
}

/**
 * @brief Set bInt vector
 *        The inner vector of b, where AInt * xInt = bInt
 *
 * @param   xVec    	x in vector form
 * @param   bVec    	b in vector form
 */
void PoissonSolver::SetVectors(double* xVec, double* bVec) {

    // Get interior matrix from bVec
    for (int i = 0; i < bIntNx; i++){
        for (int j = 0; j < bIntNy; j++) {
            bInt[i*bIntNy+j] = bVec[(i+1)*Ny+(j+1)];
        }
    }
    
    // Apply boundary conditions at edges of bInt
    cblas_daxpy(bIntNy, seconddiag, xVec + 1, 1, bInt, 1);                                  // Left BC
    cblas_daxpy(bIntNy, seconddiag, xVec + (Nx-1)*Ny + 1, 1, bInt + bIntNy*(bIntNx-1), 1);  // Right BC
    cblas_daxpy(bIntNx, firstdiag, xVec + 2*Ny - 1, Ny, bInt + (bIntNy-1), bIntNy);         // Top BC
    cblas_daxpy(bIntNx, firstdiag, xVec + Ny, Ny, bInt, bIntNy);                    	    // Bottom BC
}

/**
 * @brief Initialize BLACS for ScaLaPack
 *
 * @param   px    Number of partitions in x-direction 
 * @param   py    Number of partitions in y-direction 
 */
void PoissonSolver::InitialiseScaLaPack(int px, int py) {

    Px = px;
    Py = py;
    int MPI_init;

    // Check if MPI was initialised
    MPI_Initialized(&MPI_init);
    if (!MPI_init) {
        cout << "Error: MPI not initialised" << endl;
        throw exception();
    } else {
        // Get comm rank and size of each process
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    }
    
    char order = 'C';  // block cyclic column-major processor mapping
    // Initialise BLACS world communicator (calls MPI_Init if needed)
    Cblacs_pinfo(&mype, &npe);
    // Gets default system context (i.e. MPI_COMM_WORLD)
    Cblacs_get(0, 0, &ctx);
    // Initialise a process grid of 1 row and Px*Py columns
    Cblacs_gridinit(&ctx, &order, 1, Px * Py);
    // Gets grid info to verify the set up is complete
    Cblacs_gridinfo(ctx, &nrow, &ncol, &myrow, &mycol);
}

/**
 * @brief Generate AInt matrix in a banded storage format for ScaLaPack
 *        Ex. of column-major storage of Laplacian banded symmetric matrix for innerNx = 3 & innerNy = 3
 *         | *  *  *  *  *  *  *  *  * |
 *         | *  *  *  *  *  *  *  *  * |   KL padding = 3
 *         | *  *  *  *  *  *  *  *  * |
 *         | *  *  *  *  *  *  *  *  * |
 *         | *  *  *  *  *  *  *  *  * |   KU padding = 3
 *         | *  *  *  *  *  *  *  *  * |
 *         |  *  *  * -s -s -s -s -s -s|
 *         | *  *  *  *  *  *  *  *  * |
 *         |-f -f  * -f -f  * -f -f  * |
 *         | d  d  d  d  d  d  d  d  d |
 *         | * -f -f  * -f -f  * -f -f |
 *         | *  *  *  *  *  *  *  *  * |
 *         |-s -s -s -s -s -s  *  *  * |    
 *         where the size of the matrix for dgbsv is (1 + 2KL + 2KU) × N, KL = KU = 3, N = innerNx * innerNy
 */
void PoissonSolver::GenerateScaLaPackAIntMatrix() {

    // Define AInt banded storage matrix dimensions
    AIntNx = bIntNx * bIntNy;
    AIntNy = bIntNy * 4 + 1;

    // Initialise zeros interior A matrix
    AInt = new double[AIntNx * AIntNy]{};
    
    for (int i = 0; i < AIntNx; i++) {
        
        // Set main diagonal band
        AInt[i*AIntNy + 3 * bIntNy] = diag;

        // Set lower first super diagonal band
        if ((i+1) % bIntNy != 0) {
            AInt[i*AIntNy + 1 + 3 * bIntNy] = -1 * firstdiag;
        }

        // Set upper first super diagonal band
        if (i % bIntNy != 0) {
            AInt[i*AIntNy + 3 * bIntNy - 1] = -1 * firstdiag;
        }

        // Set lower second super diagonal band
        if (i < bIntNy * (bIntNx - 1)) {
            AInt[i*AIntNy + 4 * bIntNy] = -1 * seconddiag;
        }

        // Set upper second super diagonal band
        if (i > bIntNy - 1) {
            AInt[i*AIntNy + 2 * bIntNy] = -1 * seconddiag;
        }
    }
}

/**
 * @brief Generate AInt matrix in a banded storage format for LaPack
 *        Ex. of column-major storage of Laplacian banded symmetric matrix for innerNx = 3 & innerNy = 3
 *         | *  *  *  *  *  *  *  *  * |
 *         | *  *  *  *  *  *  *  *  * |   KL padding = 3
 *         | *  *  *  *  *  *  *  *  * |
 *         |  *  *  * -s -s -s -s -s -s|
 *         | *  *  *  *  *  *  *  *  * |
 *         |-f -f  * -f -f  * -f -f  * |
 *         | d  d  d  d  d  d  d  d  d |
 *         | * -f -f  * -f -f  * -f -f |
 *         | *  *  *  *  *  *  *  *  * |
 *         |-s -s -s -s -s -s  *  *  * |    
 *         where the size of the matrix for dgbsv is (1 + 2KL + KU) × N, KL = KU = 3, N = innerNx * innerNy
 */
void PoissonSolver::GenerateLaPackAIntMatrix() {

    // Define AInt banded storage matrix dimensions
    AIntNx = bIntNx * bIntNy;
    AIntNy = bIntNy * 3 + 1;

    // Initialize zeros interior A matrix
    AInt = new double[AIntNx * AIntNy]{};
    
    for (int i = 0; i < AIntNx; i++) {
        
        // Set main diagonal band
        AInt[i*AIntNy + 2 * bIntNy] = diag;

        // Set lower first super diagonal band
        if ((i+1) % bIntNy != 0) {
            AInt[i*AIntNy + 1 + 2 * bIntNy] = -1 * firstdiag;
        }

        // Set upper first super diagonal band
        if (i % bIntNy != 0) {
            AInt[i*AIntNy + 2 * bIntNy - 1] = -1 * firstdiag;
        }

        // Set lower second super diagonal band
        if (i < bIntNy * (bIntNx - 1)) {
            AInt[i*AIntNy + 3 * bIntNy] = -1 * seconddiag;
        }

        // Set upper second super diagonal band
        if (i > bIntNy - 1) {
            AInt[i*AIntNy + bIntNy] = -1 * seconddiag;
        }
    }
}

/**
 * @brief Update interior values of a given matrix with xInt
 *
 * @param   xVec    Full sized matrix to be updated
 */
void PoissonSolver::Updatex(double* xVec) {

    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            xVec[i*Ny+j] = xInt[(i-1)*bIntNy + (j-1)];
        }
    }
}

/**
 * @brief Set AInt matrix and its dimensions
 *
 * @param  aint    Pointer to AInt matrix
 * @param  aintnx  Grid points in x-direction
 * @param  aintny  Grid points in y-direction
 */
void PoissonSolver::SetScaLaPackAIntMatrix(double* aint, int aintnx, int aintny) {
    AInt = aint;
    AIntNx = aintnx;
    AIntNy = aintny;
}

/**
 * @brief Pre-factor AInt matrix for parallel solve
 */
void PoissonSolver::PrefactorAIntMatrixParallel() {

    // Define variables for pdgbtrf
    n = bIntNy * bIntNx;  // (global) Total problem size
    nb = max(bIntNy * 2 + 1, (int) ceil((double) n / (double) npe)); // Blocking size (number of columns per process)
    bwl = bIntNy;         // (global) Lower bandwidth
    bwu = bIntNy;         // (global) Upper bandwidth
    ja = 1;               // (global) Start offset in matrix (fortran starts at 1)
    int info;             // Variable to store success
        
    // Fill desca
    desca[0] = 501;             // banded matrix (1-by-P process grid)
    desca[1] = ctx;             // Context
    desca[2] = n;               // (global) Size of matrix being distributed
    desca[3] = nb;              // Blocking of matrix
    desca[4] = 0;               // Process column of first column of distributed A
    desca[5] = 1+2*bwl+2*bwu;   // Local leading dim
    desca[6] = 0;               // Reserved
    
    // Fill local AInt with local partition of global AInt
    prefactoredAInt = new double[AIntNy * nb]{};
    int offseta = nb * mype * AIntNy;

    for (int i = 0; i < nb; i++) {
        for (int j = 0; j < AIntNy; j++) {
            if ((i*AIntNy + j + offseta) < (AIntNy * AIntNx)) {
                prefactoredAInt[i*AIntNy + j] = AInt[i*AIntNy + j + offseta];
            }
        }
    }

    // Generate af
    laf = (nb+bwu)*(bwl+bwu)+6*(bwl+bwu)*(bwl+2*bwu);   // Given by ScaLaPack documentation
    af = new double[laf];
    
    // Query for optimal workspace. Since lwork minimally 1, this guarantees the function call fails and returns optimal lwork. Workopt is where this value will be stored in
    int lwork = -1;
    double workopt;
    
    F77NAME(pdgbtrf)(n, bwl, bwu, prefactoredAInt, ja, desca, ipiv, af, laf, &workopt, lwork, info);
    
    // Allocate optimal workspace
    lwork = (int) workopt;
    double* work = new double[lwork];
    
    // Generate ipiv. Minimum size of nb+bwu+bwl (bug described in forum http:/icl.cs.utk.edu/lapack-forum/viewtopic.php?f=13&t=2243)
    // Original documentation gives nb, but this was found to be insufficient with larger array sizes
    ipiv = new int[nb+bwu+bwl];
    
    // Factorise matrix via LU decomposition
    F77NAME(pdgbtrf)(n, bwl, bwu, prefactoredAInt, ja, desca, ipiv, af, laf, work, lwork, info);
    
    // Check if error occurred in pre-factoring
    if (info) {
        cout << "Error occurred in PDGBTRF: " << info << endl;
    }
    
    // Deallocate temp work array
    delete[] work;
}

/**
 * @brief Solve xInt in equation AInt * xInt = bInt in parallel
 */
void PoissonSolver::SolveParallel() {
    
    // Define variables for pdgbtrs
    const int nrhs = 1; // Number of RHS to solve
    double* localbint;  // Local array b
    const int ib = 1;   // Start offset in RHS vector (fortran starts at 1)
    int descb[7];       // Descriptor for RHS
    int info;           // Status value
    
    // Fill local bInt with local partition of global bInt
    localbint = new double[nb]{};
    int offsetb = nb * mype;

    for (int i = 0; i < nb; i++) {
        if ((i + offsetb) < (bIntNx * bIntNy)) {
            localbint[i] = bInt[i + offsetb];
        }
    }
    
    // Fill descb
    descb[0] = 502; // Type
    descb[1] = ctx; // Context
    descb[2] = n;   // Problem size
    descb[3] = nb;  // Blocking of matrix
    descb[4] = 0;   // Process row/column
    descb[5] = nb;  // Local leading dim
    descb[6] = 0;   // Reserved
    
    // Allocate workspace
    int lwork = (nb+bwu)*(bwl+bwu)+6*(bwl+bwu)*(bwl+2*bwu); // Given by ScaLaPack documentation
    double* work = new double[lwork];
    
    // Solve bInt
    F77NAME(pdgbtrs)('N', n, bwl, bwu, nrhs, prefactoredAInt, ja, desca, ipiv, localbint, ib, descb, af, laf, work, lwork, info);
    
    // Checke if error occured during solve
    if (info) {
        cout << "Error occurred in PDGBTRS: " << info << endl;
    }
    
    if (mype > 0) {
        // Send calculated local bInt to root process
        int rdest = 0;
        int cdest = 0;
        F77NAME(dgesd2d)(ctx, nb, 1, localbint, nb, rdest, cdest);
    } else {
        // Receive calculated local bInt from all other processes
        double* assembledB = new double[nb * npe]{};
        cblas_dcopy(nb, localbint, 1, assembledB, 1);

        for (int src = 1; src < npe; src++) {
            F77NAME(dgerv2d)(ctx, nb, 1, assembledB + src * nb, nb, 0, src);
        }
        
        // Assemble global bInt from individual local segments
        cblas_dcopy(bIntNx * bIntNy, assembledB, 1, xInt, 1);
        delete[] assembledB;
    }
    
    // Clean up memory
    delete[] work;
    delete[] localbint;
}


/**
 * @brief Pre-factor AInt matrix for serial solve
 */
void PoissonSolver::PrefactorAIntMatrixSerial() {

    // Define variables for dgbtrf
    m = bIntNy * bIntNx;      // Number of rows in matrix A
    n = bIntNy * bIntNx;      // Number of columns in matrix A
    bwl = bIntNy;             // Lower diagonal bandwidth
    bwu = bwl;                // Upper diagonal bandwidth
    lda = 1 + 2 * bwl + bwu;  // Number of rows in compressed matrix
    int info;                 // Variable to store success

    // Generate ipiv
    ipiv = new int[n];
    
    // Factorise matrix via LU decomposition
    F77NAME(dgbtrf) (m, n, bwl, bwu, AInt, lda, ipiv, info);
    
    // Checke if error occurred when prefactoring
    if (info) {
        cout << "Error occurred in DGBTRF: " << info << endl;
    }
}

/**
 * @brief Solve xInt in equation AInt * xInt = bInt in serial
 */
void PoissonSolver::SolveSerial() {

    // Define variables for dgbtrs
    int info;
    int nrhs = 1;
    
    // Solve bInt
    cblas_dcopy(bIntNy * bIntNx, bInt, 1, xInt, 1);
    F77NAME(dgbtrs) ('N', n, bwl, bwu, nrhs, AInt, lda, ipiv, xInt, n, info);
    
    // Check if error occured during solve
    if (info) {
        cout << "Error occurred in DGBTRS: " << info << endl;
    }
}
