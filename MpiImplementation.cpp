/**
 * @file MpiImplementation.cpp
 *
 * High-Performance Computing
 *
 * MPI implementation to the Coursework Assignment
 *
 * Contain functions to control MPI functionalities, crucial to controlling processor interdependencies
 *
 */
#include "MpiImplementation.h"
#include <iostream>
#include <mpi.h>

using namespace std;

/**
 * @brief Check if number of processors is the same as total number of partitions
 *
 * @param   np     Number of processors
 * @param   Px     Number of partitions in the x-direction (parallel)
 * @param   Py     Number of partitions in the y-direction (parallel)
 */
bool Processors(const int &np, const int &Px, const int &Py) {
	if (np == Px * Py) {
		return true;
    	} else {
		return false;
    	}
}

/**
 * @brief Determine neighbours in all 4 directions for current rank
 *
 * @param   gridcomm       	New communicator based on Cartesian Topology
 * @param   neighbour_rank	Neighbours at top, bottom, left and right of a rank accordingly
 */ 
void Neighbours(MPI_Comm gridcomm, int* neighbour_rank) {
	MPI_Cart_shift(gridcomm, 0, -1, &neighbour_rank[0], &neighbour_rank[2]);	// top & bottom
	MPI_Cart_shift(gridcomm, 1,  1, &neighbour_rank[1], &neighbour_rank[3]);	// left & right
}

/**
 * @brief Distribution of workload as evenly as possible to each processor
 *
 * @param   rank         	Rank of a process within group of comm (gridcomm)
 * @param   dx           	Grid spacing in x-direction
 * @param   dy           	Grid spacing in y-direction
 * @param   cart_coords       	Array containing the Cartesian coordinates of specified process
 * @param   global_coords       Array containing the global coordinates of specified process
 * @param   nx           	Number of grid points of one partition in x-direction
 * @param   ny           	Number of grid points of one partition in y-direction
 */  
void Allocation(const int &rank, const int &Nx, const int &Ny, const int &Px, const int &Py, const double &dx, const double &dy, int* cart_coords, double* global_coords, int &nx, int &ny) {
    
    // Ranks with smaller coordinates have 1 larger n_dim in grids
    nx = (Nx - 2) / Px;
    ny = (Ny - 2) / Py;
    int remainder_x = (Nx - 2) % Px;
    int remainder_y = (Ny - 2) % Py;

    // global_coords[i, j] and cart_coords[j, i]. This is due to the axis defined in the problem.  
    // x direction
    if (cart_coords[1] < remainder_x) {
        nx++;
        global_coords[0] = (cart_coords[1] * nx + 1)*dx;
    } else {
	global_coords[0] = (cart_coords[1] * nx + 1)*dx;
    }

    // y direction
    if (cart_coords[0] < remainder_y) {
        ny++;
        global_coords[1] = (cart_coords[0] * ny + 1) * dy;
    } else {
	global_coords[1] = (cart_coords[0] * ny + 1) * dy;
    }
}