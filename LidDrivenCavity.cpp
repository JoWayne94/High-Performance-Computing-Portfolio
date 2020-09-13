/**
 * @file LidDrivenCavity.cpp
 *
 * High-Performance Computing
 *
 * Solution to the Coursework Assignment
 *
 * Solves the lid-driven cavity problem using the algorithm
 */
#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include "cblas.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <mpi.h>

using namespace std;

/**
 * @brief Constructor
 *
 * @param   gridcomm       	New communicator based on Cartesian Topology
 * @param   rank         	Rank of a process within group of communicator (gridcomm)
 * @param   cart_coords       	Array containing the Cartesian coordinates of specified process
 * @param   global_coords	Array containing the global coordinates of specified process
 * @param   neighbour_rank	Neighbours at top, bottom, left and right of a rank accordingly
 * @param   dt_bool	 	Checker of the restriction on the time-step
 * @param   deltat       	Time step size
 * @param   finalt       	Final time
 * @param   nx        	 	Number of grid points in x-direction
 * @param   ny        	 	Number of grid points in y-direction
 * @param   re           	Reynolds number
 */
LidDrivenCavity::LidDrivenCavity(MPI_Comm gridcomm, int rank, int *cart_coords, double *global_coords, int *neighbour_rank, bool &dt_bool, double deltat, double finalt, int nx, int ny, double re, double dx, double dy) {

    // Set grid size, final time, reynolds number, and time step
    dt_bool = 1;
    T = finalt;
    Nx = nx;		// Initialize Nx and Ny arrays
    Ny = ny;
    Re = re;

    n_dim = Nx * Ny;    // number of elements in a single partition
    ld_A = Ny + 1;      // leading dimension of the A symmetric banded matrix
    ld_B = 3;	        // leading dimension of the y denominator banded matrix at LHS of interior vorticity at time t + dt step
    ld_C = Ny * 2 + 1;  // leading dimension of the x denominator banded matrix at LHS of interior vorticity at time t + dt step

    this-> gridcomm = gridcomm;
    this-> rank = rank;
    this-> dx = dx;
    this-> dy = dy;

    // Copy memory blocks
    memcpy(this-> cart_coords, cart_coords, 2 * sizeof(int));
    memcpy(this-> global_coords, global_coords, 2 * sizeof(double));
    memcpy(this-> neighbour_rank, neighbour_rank, 4 * sizeof(int));

   // Restriction on time-step is given by dt < Re*dx*dy/4
   try
   {
	dt = deltat;
	if (dt >= Re*dx*dy/4) {
		throw std::logic_error("Courant-Friedrichs-Lewy condition not satisfied, try new inputs.");
	}
    }
   catch (const std::logic_error &e)
    {
	if (rank == 0) {
    		cout << "An error occured: " << e.what() << endl;
		dt_bool = 0;
    	}
    }
}

/**
 * @brief Destructor
 *
 * 	  Clean up memory
 */
LidDrivenCavity::~LidDrivenCavity() {

	delete[] v;
	delete[] s;
	delete[] vbc_top;
	delete[] vbc_bottom;
	delete[] vbc_left;
	delete[] vbc_right;
	delete[] sbc_top;
	delete[] sbc_bottom;
	delete[] sbc_left;
	delete[] sbc_right;

	delete[] A;
	delete[] B;
	delete[] C;
}

/**
 * @brief Set domain size
 *	  Initialize Lx and Ly arrays
 *
 * @param   xlen        Length of the domain in the x-direction
 * @param   ylen        Length of the domain in the y-direction
 */
void LidDrivenCavity::SetDomainSize(double xlen, double ylen) {
    Lx = xlen;
    Ly = ylen;
}

/**
 * @brief Initialize all variables
 */
void LidDrivenCavity::Initialise() {

    // Vorticity and stream function matrices 
    v = new double[n_dim];
    s = new double[n_dim];
    s_prev = new double[n_dim];  // Stream function matrix of previous time-step used to calculate the error difference

    // Vorticity and stream function boundary arrays for use of boundary conditions
    vbc_top = new double[Nx];
    vbc_bottom = new double[Nx];
    vbc_left = new double[Ny];
    vbc_right = new double[Ny];
    sbc_top = new double[Nx];
    sbc_bottom = new double[Nx];
    sbc_left = new double[Ny];
    sbc_right = new double[Ny];

    // Horizontal and vertical velocities
    Ux = new double[n_dim];
    Uy = new double[n_dim];

    // Fill variables to a full zero matrix
    cblas_dscal(n_dim, 0.0, v, 1);
    cblas_dscal(n_dim, 0.0, s, 1);
    cblas_dscal(Nx, 0.0, vbc_top, 1);
    cblas_dscal(Nx, 0.0, vbc_bottom, 1);
    cblas_dscal(Ny, 0.0, vbc_left, 1);
    cblas_dscal(Ny, 0.0, vbc_right, 1);
    cblas_dscal(Nx, 0.0, sbc_top, 1);
    cblas_dscal(Nx, 0.0, sbc_bottom, 1);
    cblas_dscal(Ny, 0.0, sbc_left, 1);
    cblas_dscal(Ny, 0.0, sbc_right, 1);
}

/**
 * @brief Create banded matrices
 */
void LidDrivenCavity::CreateBandedMatrices() {

	// Banded matrices used for algorithm
	A = new double[n_dim * ld_A];
	B = new double[n_dim * ld_B];
	C = new double[n_dim * ld_C];

	// Diagonal values in the A original matrix
	double diag = 2/(dx*dx) + 2/(dy*dy);	// main diagonal values
	double first_diag = -1/(dy*dy);		// diagonal values right above and below the main diagonal
	double second_diag = -1/(dx*dx);	// diagonal values furthest away from the main diagonal

	// Create A, B, and C matrices
	for (unsigned int i = 0; i < n_dim; i++) {
		A[ld_A - 1 + i * ld_A] 			  = diag;
		if (i % Ny != 0)   A[ld_A - 2 + i * ld_A] = first_diag;
		if (i >= Ny)       A[i * ld_A] 		  = second_diag;

		if (i % Ny != 0)     B[i*ld_B] 	   = 1/(2*dy);
		if (i % Ny != Ny-1)  B[2 + i*ld_B] = -1/(2*dy);

		if (i >= Ny)       C[i*ld_C] 	      = 1/(2*dx);
		if (i < n_dim-Ny)  C[ld_C-1 + i*ld_C] = -1/(2*dx);
	}
}

/**
 * @brief Output boundary array b in the linear system y = Ax + b for matrix A
 * 
 * @param   b       	Boundary array for partitions
 * @param   x       	Input to decide either vorticity or stream function
 */
void LidDrivenCavity::BoundariesPartitionA(double *b, char x) {
	static double first_diag = -1/(dy*dy);
	static double second_diag = -1/(dx*dx);

	if (x == 's') {
		for (unsigned int i = 0; i < Ny; i++)
		{
			b[i] += sbc_left[i] * second_diag;
			b[i+(Nx-1)*Ny] += sbc_right[i] * second_diag;
		}
		for (unsigned int i = 0; i < Nx; i++)
		{
			b[i*Ny] += sbc_bottom[i] * first_diag;
			b[i*Ny + Ny-1] += sbc_top[i] * first_diag;
		}
	} else if (x == 'v') {
		for (unsigned int i = 0; i < Ny; i++)
		{
			b[i] += vbc_left[i] * second_diag;
			b[i+(Nx-1)*Ny] += vbc_right[i] * second_diag;
		}
		for (unsigned int i = 0; i < Nx; i++)
		{
			b[i*Ny] += vbc_bottom[i] * first_diag;
			b[i*Ny + Ny-1] += vbc_top[i] * first_diag;
		}
	}
}

/**
 * @brief Output boundary array b in the linear system y = Bx + b for matrix B
 * 
 * @param   b       	Boundary array for partitions
 * @param   x       	Input to decide either vorticity or stream function
 */
void LidDrivenCavity::BoundariesPartitionB(double *b, char x) {

	if (x == 's') {
		for (unsigned int i = 0; i < Nx; i++)
		{
			b[i*Ny] = sbc_bottom[i] * (-1/(2*dy));
			b[i*Ny + Ny-1] = sbc_top[i] * 1/(2*dy);
		}
	} else if (x == 'v') {
		for (unsigned int i = 0; i < Nx; i++)
		{
			b[i*Ny] = vbc_bottom[i] * (-1/(2*dy));
			b[i*Ny + Ny-1] = vbc_top[i] * 1/(2*dy);
		}
	}
}

/**
 * @brief Output boundary array b in the linear system y = Cx + b for matrix C
 * 
 * @param   b       	Boundary array for partitions
 * @param   x       	Input to decide either vorticity or stream function
 */
void LidDrivenCavity::BoundariesPartitionC(double *b, char x) {

	if (x == 's') {
		for (unsigned int i = 0; i < Ny; i++)
		{
			b[i] = sbc_left[i] * (-1/(2*dx));
			b[i+(Nx-1)*Ny] = sbc_right[i] * 1/(2*dx);
		}
	} else if (x == 'v') {
		for (unsigned int i = 0; i < Ny; i++)
		{
			b[i] = vbc_left[i] * (-1/(2*dx));
			b[i+(Nx-1)*Ny] = vbc_right[i] * 1/(2*dx);
		}
	}
}

/**
 * @brief Calculate horizontal and vertical velocity from vorticity or stream function (equation 3)
 * 	  Ux = {s_{i, j+1}^{n} - s_{i, j-1}^{n}} / {2 * dy}
 * 	  Uy = {s_{i+1, j}^{n} - s_{i-1, j}^{n}} / {2 * dx} 
 * 	  Through the Finite Difference Method
 */
void LidDrivenCavity::CalculateVelocity() {
	BoundariesPartitionB(Ux, 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, n_dim, n_dim, 1,  1,   1.0, B, ld_B, s, 1,  1.0, Ux, 1);
	
	BoundariesPartitionC(Uy, 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, n_dim, n_dim, Ny, Ny, -1.0, C, ld_C, s, 1, -1.0, Uy, 1);
}

/**
 * @brief Calculation of vorticity boundary conditions (equations 6 to 9)
 */
void LidDrivenCavity::VorticityBCs() {

    // Top boundary, equation 6
    if (neighbour_rank[0] == MPI_PROC_NULL)		// if it is -2
    {
        for (unsigned int i = 0; i < Nx; i++)
        {
            vbc_top[i] = (sbc_top[i] - s[i * Ny + (Ny - 1)]) * 2/(dy*dy) - 2*U/dy;
        }
    }
    // Bottom boundary, equation 7
    if (neighbour_rank[2] == MPI_PROC_NULL)
    {
        for (unsigned int i = 0; i < Nx; i++)
        {
            vbc_bottom[i] = (sbc_bottom[i] - s[i * Ny]) * 2/(dy*dy);
        }
    }
    // Left boundary, equation 8
    if (neighbour_rank[1] == MPI_PROC_NULL)
    {
        for (unsigned int i = 0; i < Ny; i++)
        {
            vbc_left[i] = (sbc_left[i] - s[i]) * 2/(dx*dx);
        }
    }
    // Right boundary, equation 9
    if (neighbour_rank[3] == MPI_PROC_NULL)
    {
        for (unsigned int i = 0; i < Ny; i++)
        {
            vbc_right[i] = (sbc_right[i] - s[(Nx - 1) * Ny + i]) * 2/(dx*dx);
        }
    }
}

/**
 * @brief Calculation of interior vorticity at time t (equation 10)
 */
void LidDrivenCavity::InteriorVorticity() {
    	double *b =  new double[n_dim];
	BoundariesPartitionA(b, 's');

	// Solves b = As + b. Output is b.
	cblas_dsbmv (CblasColMajor, CblasUpper, n_dim, Ny, 1.0, A, ld_A, s, 1, 1.0, b, 1);
	
	// Puts b back into vorticity, v = b. Needs to output v.
	cblas_dcopy (n_dim, b, 1, v, 1);

	// Clean up
	delete[] b;
}

/**
 * @brief Calculation of interior vorticity at time t + dt (equation 11)
 *
 * 	  v = v + (dt/Re) * b_A + dt * niv_lhs1 - dt * niv_lhs2
 */
void LidDrivenCavity::NextInteriorVorticity() {

	double* b_A = new double[n_dim];	// RHS term of equation 11
	double* niv_lhs1 = new double[n_dim];	// First LHS term of equation 11
	double* niv_lhs2 = new double[n_dim];	// Second LHS term of equation 11

    	// Fill b_A, niv_lhs1, and niv_lhs2 to a full zero matrix
	cblas_dscal(n_dim, 0.0, b_A, 1);
	cblas_dscal(n_dim, 0.0, niv_lhs1, 1);
	cblas_dscal(n_dim, 0.0, niv_lhs2, 1);

	// Solves b_A = Av + b_A
	BoundariesPartitionA(b_A, 'v');

	// Solves b = Av + b, output is b_A
	cblas_dsbmv (CblasColMajor, CblasUpper, n_dim, Ny, 1.0, A, ld_A, v, 1, 1.0, b_A, 1);

	// Calculate niv_lhs1
	double* temp_x1 = new double[n_dim];	// temp_x1 = Cs + b, b refers to the boundary term
	double* temp_y1 = new double[n_dim];	// temp_y1 = Bv + b

	cblas_dscal(n_dim, 0.0, temp_x1, 1);
	cblas_dscal(n_dim, 0.0, temp_y1, 1);

	BoundariesPartitionC(temp_x1, 's');
	BoundariesPartitionB(temp_y1, 'v');

	cblas_dgbmv (CblasColMajor, CblasNoTrans, n_dim, n_dim, Ny, Ny, 1.0, C, ld_C, s, 1, 1.0, temp_x1, 1);
	cblas_dgbmv (CblasColMajor, CblasNoTrans, n_dim, n_dim, 1,  1,  1.0, B, ld_B, v, 1, 1.0, temp_y1, 1);

	for (unsigned int i = 0; i < n_dim; i++) {
		niv_lhs1[i] = temp_x1[i] * temp_y1[i];
	}

	// Clean up
	delete[] temp_x1;
	delete[] temp_y1;

	// Calculate niv_lhs2
	double* temp_x2 = new double[n_dim];	// temp_x2 = Cv + b, b is the boundary term
	double* temp_y2 = new double[n_dim];	// temp_y2 = Bs + b

	cblas_dscal(n_dim, 0.0, temp_x2, 1);
	cblas_dscal(n_dim, 0.0, temp_y2, 1);

	BoundariesPartitionC(temp_x2, 'v');
	BoundariesPartitionB(temp_y2, 's');

	cblas_dgbmv (CblasColMajor, CblasNoTrans, n_dim, n_dim, Ny, Ny, 1.0, C, ld_C, v, 1, 1.0, temp_x2, 1);
	cblas_dgbmv (CblasColMajor, CblasNoTrans, n_dim, n_dim, 1,  1,  1.0, B, ld_B, s, 1, 1.0, temp_y2, 1);

	for (unsigned int i = 0; i < n_dim; i++) {
		niv_lhs2[i] = temp_x2[i] * temp_y2[i];
	}

	// Clean up
	delete[] temp_x2;
	delete[] temp_y2;

	// Calculate equation 11
	cblas_daxpy (n_dim, -dt/Re, b_A,      1, v, 1);
	cblas_daxpy (n_dim,  dt,    niv_lhs1, 1, v, 1);
	cblas_daxpy (n_dim, -dt,    niv_lhs2, 1, v, 1);

	// Clean up
	delete[] b_A;
	delete[] niv_lhs1;
	delete[] niv_lhs2;
}

/**
 * @brief Solution for the Poisson problem to compute the stream-function at time t + dt (equation 12)
 * 
 * 	  Create a new PoissonSolver instance and solve Poisson
 */
void LidDrivenCavity::Integrate() {

	if (Poisson_Solver == nullptr)
	{
		Poisson_Solver = new PoissonSolver(Nx, Ny, dx, dy);
	}
	Poisson_Solver-> BoundariesPartition(sbc_top, sbc_bottom, sbc_left, sbc_right);
	Poisson_Solver-> SolvePoisson(s, v);
}

/**
 * @brief Send and receive inputs and outputs between processors
 * 
 * @param   x        	Either vorticity or stream function
 * @param   x_top    	Top boundary array of a partition
 * @param   x_bottom    Bottom boundary array of a partition
 * @param   x_left   	Left most boundary array of a partition
 * @param   x_right  	Right most boundary array of a partition
 */
void LidDrivenCavity::SendRecvProcess(double* x, double* x_top, double* x_bottom, double* x_left, double* x_right) {
    
    // Left partition (right most array) to Right partition (left most array)
    // Right partition (left most array) to Left partiton (right most array)
    MPI_Sendrecv(x + (Nx-1)*Ny, Ny, MPI_DOUBLE, neighbour_rank[3], 3, x_left,  Ny, MPI_DOUBLE, neighbour_rank[1], 3, gridcomm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(x, 		Ny, MPI_DOUBLE, neighbour_rank[1], 1, x_right, Ny, MPI_DOUBLE, neighbour_rank[3], 1, gridcomm, MPI_STATUS_IGNORE);

    // Top to bottom and bottom to top
    // Pointers cannot be multiplied, therefore temporary arrays are needed
    // Initialize temporary bottom and top boundary arrays
    double bottom[Nx];
    double top[Nx];

    cblas_dcopy(Nx, x, Ny, bottom, 1);
    MPI_Sendrecv(bottom, Nx, MPI_DOUBLE, neighbour_rank[2], 2, x_top,    Nx, MPI_DOUBLE, neighbour_rank[0], 2, gridcomm, MPI_STATUS_IGNORE);

    cblas_dcopy(Nx, x + Ny-1, Ny, top, 1);
    MPI_Sendrecv(top,    Nx, MPI_DOUBLE, neighbour_rank[0], 0, x_bottom, Nx, MPI_DOUBLE, neighbour_rank[2], 0, gridcomm, MPI_STATUS_IGNORE);
}

/**
 * @brief Run through entire algorithm over the specified time period, t = [0, T].
 *
 * 	  Stop either at the end of time domain or when a steady state solution is computed.
 */
void LidDrivenCavity::Algorithm(){

    // Initialize termination criteria
    double t = 0.0;
    double tol = 0.0001; 	// tolerance of the maximum residual
    double error = 0.0;		// residual, take 2-norm of stream function
    double error_max = 0.0;	// maximum residual along the whole domain
    double s_norm = 0.0;	// 2-norm of stream function

    // Configure the solver here
    Initialise();
    SetDomainSize(Lx, Ly);
    CreateBandedMatrices();

    while (t <= T) {

        cblas_dcopy(n_dim, s, 1, s_prev, 1);
        VorticityBCs();
        InteriorVorticity();
        SendRecvProcess(v, vbc_top, vbc_bottom, vbc_left, vbc_right);
        NextInteriorVorticity();

	// 5 stream function arrays, therefore iterate 5 times
        for (unsigned int i = 0; i < 5; i++) {
            Integrate();
            SendRecvProcess(s, sbc_top, sbc_bottom, sbc_left, sbc_right);
        }

		t += dt;
        s_norm = cblas_dnrm2(n_dim, s, 1);
        cblas_daxpy(n_dim, -1.0, s, 1, s_prev, 1);
        error = cblas_dnrm2(n_dim, s_prev, 1);  	// calculation of 2-norm with cblas routine
        error /= s_norm;

	// Find maximum residual in all errors, outputs max value to error_max
        MPI_Reduce(&error, &error_max, 1, MPI_DOUBLE, MPI_MAX, 0, gridcomm);
        if (rank == 0) {
                cout << "t = " << setw(5) << t << setw(30) << "Maximum residual is " << error_max << endl;
        }

	// Broadcast maximum error to all processes
        MPI_Bcast(&error_max, 1, MPI_DOUBLE, 0, gridcomm);
        if(error_max < tol) {
		if (rank == 0) {
			cout << "Termination criteria satisfied, exit program." << endl;
		} break;
	}
    }

    // Finally, calculate the velocities.
    CalculateVelocity();
}

#define FORMAT(vOut, x, y, v, s, Ux, Uy) vOut << setw(10) << x << setw(10) << y << setw(12) << v << setw(12) << s << setw(12) << Ux << setw(12) << Uy << endl;

// Outputs to Output.txt, all data to be used to plot on Matlab.
void LidDrivenCavity::Output(int Px, int Py, double Lx, double Ly)
{
	for (int k = 0; k < Px *Py; k++)
	{
		if (k == rank)
		{
			// Write header and corner values. Corner values are known to be 0 all the time.
			if (k == 0)
			{
				ofstream vOut("Output.txt", ios::out | ios::trunc);
				FORMAT(vOut, "x", "y", "v", "s", "Ux", "Uy");
				FORMAT(vOut, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
				FORMAT(vOut, 0.0, Ly, 0.0, 0.0, 0.0, U);
				FORMAT(vOut, Lx, 0.0, 0.0, 0.0, 0.0, 0.0);
				FORMAT(vOut, Lx, Ly, 0.0, 0.0, 0.0, 0.0);
				vOut.close();
			}
			ofstream vOut("Output.txt", ios::out | ios::app);
			double x, y;
			// Write interior values of v, s, Ux and Uy
			vOut.precision(5);
			for (int i = 0; i < Nx; i++)
			{
				for (int j = 0; j < Ny; j++)
				{
					x = global_coords[0] + i*dx;
					y = global_coords[1] + j*dy;
					FORMAT(vOut, x, y, v[i*Ny + j], s[i*Ny + j],
						       	Ux[i*Ny + j], Uy[i*Ny + j]);
				}
			}
			// Write boundary values

			// top. Non zero Ux value!
			if (cart_coords[0] == Py - 1)
			{
				y = Ly;
				for (int i = 0; i < Nx; i++)
				{
					x = global_coords[0] + i*dx;
					FORMAT(vOut, x, y, vbc_top[i], 0.0, U, 0.0);
				}
			}
			// bottom
			if (cart_coords[0] == 0)
			{
				y = 0;
				for (int i = 0; i < Nx; i++)
				{
					x = global_coords[0] + i*dx;
					FORMAT(vOut, x, y, vbc_bottom[i], 0.0, 0.0, 0.0);
				}
			}
			// left
			if (cart_coords[1] == 0)
			{
				x = 0;
				for (int j = 0; j < Ny; j++)
				{
					y = global_coords[1] + j*dy;
					FORMAT(vOut, x, y, vbc_left[j], 0.0, 0.0, 0.0);
				}
			}

			// right
			if (cart_coords[1] == Px - 1)
			{
				x = Lx;
				for (int j = 0; j < Ny; j++)
				{
					y = global_coords[1] + j*dy;
					FORMAT(vOut, x, y, vbc_right[j], 0.0, 0.0, 0.0);
				}
			}
			vOut.close();
		}
		MPI_Barrier(gridcomm);
	}
}