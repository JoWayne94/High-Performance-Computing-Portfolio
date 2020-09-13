/**
 * @file LidDrivenCavitySolver.cpp
 *
 * High-Performance Computing
 *
 * Solution to the Coursework Assignment
 *
 * Main file of the solution
 *
 */
#include "LidDrivenCavity.h"
#include "ProgramInputs.h"
#include "MpiImplementation.h"

#include <iostream>
#include <exception>
#include <iomanip>
#include <cmath>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv) {

    int rank;
    int n_dim;

    // Initialise MPI
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        cout << "Initialization of MPI failed." << endl;
        return 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 	// Get rank
    MPI_Comm_size(MPI_COMM_WORLD, &n_dim); 	// Get number of processors

    // Parameters from command input. Variables have been parsed through ProgramInputs.cpp
    double dt;
    double T;
    int Nx;
    int Ny;
    double Lx; 
    double Ly;
    double Re;
    int Px;
    int Py;

    po::variables_map vm;
    // Input status, 1 for successful input.
    bool status;
    status = InputStatus(argc, argv, vm);

    // Terminate the program if error occurs or help is called
    if (!status) {
        MPI_Finalize();
        return 0;
    }

    Inputs(vm, dt, T, Nx, Ny, Lx, Ly, Re, Px, Py);
    double dx = Lx/(Nx - 1);
    double dy = Ly/(Ny - 1);

    // Check if number of processors satisfy the condition, to avoid erronenous outputs.
    if (!Processors(n_dim, Px, Py)) {
        if (rank == 0) {
            cout << "Number of processors, np, does not match Px * Py. Please try again." << endl;
        }
        MPI_Finalize();
        return 0;
    }

    // Pre-operation
    MPI_Comm gridcomm;
    const int dims = 2;
    int n_dims[dims] = {Py, Px};
    int period[dims] = {0, 0};
    int reorder = 0;
    
    // New communicator based on Cartesian Topology
    MPI_Cart_create(MPI_COMM_WORLD, dims, n_dims, period, reorder, &gridcomm);
    int cart_coords[dims];

    // Assign a rank to each coordinate
    MPI_Cart_coords(gridcomm, rank, dims, cart_coords);
    int neighbour_rank[4] = {0};

    // Obtain neighborhood information for each rank
    Neighbours(gridcomm, neighbour_rank);
    int nx;				// number of grids for each rank
    int ny;	
    double global_coords[2] = {0.0, 0.0};

    // Distribution of work to each process, calling function from MpiImplementation.cpp
    Allocation(rank, Nx, Ny, Px, Py, dx, dy, cart_coords, global_coords, nx, ny);
    bool dt_bool;

    // Create a new instance of the LidDrivenCavity class for every rank
    LidDrivenCavity* solver = new LidDrivenCavity(gridcomm, rank, cart_coords, global_coords, neighbour_rank, dt_bool, dt, T, nx, ny, Re, dx, dy);
    
    // Condition in order to enforce stability of solution 
    if (!dt_bool) {
        MPI_Finalize();
        return 0;
    }

    // Run the solver
    solver-> Algorithm();
    cout << "Solver completed." << endl;

    // Output the result. Write to output.txt file
    solver-> Output(Px, Py, Lx, Ly);

    MPI_Finalize();
    return 0;
}




