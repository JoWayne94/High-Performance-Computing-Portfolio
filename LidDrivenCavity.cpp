/**
 * @file LidDrivenCavity.cpp
 *
 * High-Performance Computing
 *
 * Solution to the Coursework Assignment
 *
 * Solves the lid-driven cavity problem using various algorithms
 * 
 */
#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include "cblas.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>

using namespace std;

/**
 * @brief Constructor for the LidDrivenCavity class
 */
LidDrivenCavity::LidDrivenCavity() {
}

/**
 * @brief Destructor for the LidDrivenCavity class
 *
 * 	  Clean up memory
 */
LidDrivenCavity::~LidDrivenCavity() {

	delete[] v;
	delete[] s;
}

/**
 * @brief Initialise and configure MPI, segment the v and s array (in order of storage position) into evenly distributed chunks for each process
 */
void LidDrivenCavity::SetMPIConfig() {

    // Check if MPI was initialised
    int MPI_init;
    MPI_Initialized(&MPI_init);

    if (!MPI_init) {
        cout << "Error: MPI not initialised" << endl;
        throw exception();
    } else {

        // Get comm rank and size of each process
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        
        // Check if running serial or parallel
        if (nprocs > 1) {

            // Calculate number of elements to allocate for each partition
            int spacing = (int) floor((double) innerNxy / (double) nprocs);
            int remainder = (int)innerNxy % nprocs;
            int tag = 0;
            int source = 0;
            
            // Generate coordinate pairs and assign to partitions from root
            if (rank == source) {
                // Generate i,j interior coordinate pairs
                int* icoord = new int[innerNxy];
                int* jcoord = new int[innerNxy];

                for (unsigned int i = 0; i < innerNx; i++) {
                    for (unsigned int j = 0; j < innerNy; j++) {
                        icoord[i*innerNy+j] = i + 1;
                        jcoord[i*innerNy+j] = j + 1;
                    }
                }
                
                // Calculate start and end points of the partitions
                int* startPoints = new int [nprocs];
                int* endPoints = new int [nprocs];
                
                int currentLoc = 0;  // Current location
                for (int i = 0; i < nprocs; i++) {
                    startPoints[i] = currentLoc;
                    currentLoc += spacing;

                    if (remainder > 0) {
                        currentLoc++;
                        remainder--;
                    }
                    endPoints[i] = currentLoc - 1;
                }
                
                // Send array of partitioned i,j coordinates to different processes
                for (int destination = 0; destination < nprocs; destination++) {

                    int start = startPoints[destination];
                    int end = endPoints[destination];
                    int arrlength = end - start + 1;

                    // Split inner coordinates of interior matrix into respective partitions
                    if (destination == 0) {
                        arrlenCoord = arrlength;
                        icoordInner = new int[arrlenCoord];
                        jcoordInner = new int[arrlenCoord];
                        copy(icoord + start, icoord + end + 1, icoordInner);
                        copy(jcoord + start, jcoord + end + 1, jcoordInner);
                    } else {
                        // Create array containing the local partition of coordinates
                        int* tempicoord = new int[arrlength];
                        int* tempjcoord = new int[arrlength];
                        copy(icoord + start, icoord + end + 1, tempicoord);
                        copy(jcoord + start, jcoord + end + 1, tempjcoord);
                        
                        // Send array containing coordinates to individual processes
                        MPI_Send(&arrlength, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
                        MPI_Send(tempicoord, arrlength, MPI_INT, destination, tag, MPI_COMM_WORLD);
                        MPI_Send(tempjcoord, arrlength, MPI_INT, destination, tag, MPI_COMM_WORLD);

			// Deallocate temp arrays
                        delete[] tempicoord;
                        delete[] tempjcoord;
                    }
                }
                
                // Clean up memory
                delete[] startPoints;
                delete[] endPoints;
                delete[] icoord;
                delete[] jcoord;

            } else {
                // Receive partitioned i,j coordinates from root
                MPI_Recv(&arrlenCoord, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                icoordInner = new int[arrlenCoord];
                jcoordInner = new int[arrlenCoord];
                MPI_Recv(icoordInner, arrlenCoord, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(jcoordInner, arrlenCoord, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
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
 * @brief Set grid size
 *	  Initialize Nx and Ny arrays, and other variables
 *
 * @param   nx        Number of grid points in the x-direction
 * @param   ny        Number of grid points in the y-direction
 */
void LidDrivenCavity::SetGridSize(unsigned int nx, unsigned int ny) {
    Nx = nx;
    Ny = ny;
    
    // Matrix size of entire domain
    Nxy = Nx * Ny;
    
    // Interior matrix dimensions
    innerNx = Nx - 2;
    innerNy = Ny - 2;
    
    // Interior matrix size
    innerNxy = innerNx * innerNy;
}

/**
 * @brief Set partition size
 *	  Initialize Px and Py arrays
 *
 * @param   px        Number of partitions in the x-direction
 * @param   py        Number of partitions in the y-direction
 */
void LidDrivenCavity::SetPartitionSize(unsigned int px, unsigned int py) {
    Px = px;
    Py = py;
}

/**
 * @brief Set time step
 *
 * @param   deltat        Time step size for each iteration
 */
void LidDrivenCavity::SetTimeStepSize(double deltat) {
    dt = deltat;
}

/**
 * @brief Set final time
 *
 * @param   finalt        Final time
 */
void LidDrivenCavity::SetFinalTime(double finalt) {
    T = finalt;
}

/**
 * @brief Set Reynolds number
 *
 * @param   re        Reynolds number
 */
void LidDrivenCavity::SetReynoldsNumber(double re) {
    Re = re;
}

/**
 * @brief Set grid spacing
 *
 * @param   deltax        Grid spacing in the x-direction
 * @param   deltay        Grid spacing in the y-direction
 */
void LidDrivenCavity::SetGridSpacing(double deltax, double deltay) {
    dx = deltax;
    dy = deltay;
}

/**
 * @brief Set vorticity and stream function
 *	  Initialize v and s matrices as zeros arrays to satisfy initial conditions
 */
void LidDrivenCavity::Initialise() {
    // t = 0, omega(x,y) = 0, psi(x,y) = 0
    v = new double[Nxy]{};
    s = new double[Nxy]{};
}

/**
 * @brief Calculation of vorticity boundary conditions (equations 6 to 9)
 */
void LidDrivenCavity::SetVorticityBCs() {

    if(rank == 0) {
        for (unsigned int i = 0; i < Nx; i++) {
            // Top boundary condition, equation 6
            v[i*Ny+(Ny-1)] = (s[i*Ny+(Ny-1)] - s[i*Ny+(Ny-2)]) * 2/dy/dy - 2*U/dy;

            // Bottom boundary condition, equation 7
            v[i*Ny] = (s[i*Ny] - s[i*Ny+1]) * 2/dy/dy;
        }

        for (unsigned int i = 0; i < Ny; i++) {
            // Left boundary condition, equation 8
            v[i] = (s[i] - s[i+Ny]) * 2/dx/dx;

            // Right boundary condition, equation 9
            v[i+Ny*(Nx-1)] = (s[i+Ny*(Nx-1)] - s[i+Ny*(Nx-2)]) * 2/dx/dx;
        }
    }
}

/**
 * @brief Calculation of interior vorticity at time t (equation 10)
 *        Second order central difference scheme for the spatial derivatives
 */
void LidDrivenCavity::SetInteriorVorticity() {

    // Parallel code
    if (nprocs > 1) {
        // Initialise zeros temp matrix to store local partition of matrix
        double* vTemp = new double[Nxy]{};
        
        // Calculate interior vorticity of the assigned segment on the given process
        for (int k = 0; k < arrlenCoord; k++) {
            int i = icoordInner[k];
            int j = jcoordInner[k];
            vTemp[i*Ny+j] = -(s[i*Ny+(j+1)] - 2*s[i*Ny+j] + s[i*Ny+(j-1)])/dy/dy - (s[(i+1)*Ny+j] - 2*s[i*Ny+j] + s[(i-1)*Ny+j])/dx/dx;
        }

        if (rank > 0) {
            // Send calculated local v matrix to root process
            MPI_Send(vTemp, Nxy, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else {
            // Receive each individual local segment and combine them on root for global
            for (int src = 1; src < nprocs; src++) {
                double* vTemploc = new double[Nxy]{};
                MPI_Recv(vTemploc, Nxy, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                cblas_daxpy(Nxy, 1.0, vTemploc, 1, vTemp, 1);

		// Clean up
                delete[] vTemploc;
            }
        }
        // Update vorticity matrix
        delete[] v;
        v = vTemp;

    } else {
        // Serial code, calculate interior vorticity for entire vorticity matrix
        for (unsigned int i = 1; i < (Nx-1); i++) {
            for (unsigned int j = 1; j < (Ny-1); j++) {
                v[i*Ny+j] = -(s[i*Ny+(j+1)] - 2*s[i*Ny+j] + s[i*Ny+(j-1)])/dy/dy - (s[(i+1)*Ny+j] - 2*s[i*Ny+j] + s[(i-1)*Ny+j])/dx/dx;
            }
        }
    }   
}

/**
 * @brief Calculation of interior vorticity at time t + dt (equation 11)
 *	  Second order central difference scheme for the spatial derivatives
 */
void LidDrivenCavity::NextInteriorVorticity() {

    // Initialise matrix to store updated v. Vorticity is not updated immediately as the neighbouring segments rely on the t value to calculate t + dt
    double* temp = new double[Nxy]{};
    
    // Parallel code
    if (nprocs > 1) {
        // Broadcast original v matrix at time t to all processes
        MPI_Bcast(v, Nxy, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        // Calculate updated interior vorticity of the assigned segment on the given process
        for (int k = 0; k < arrlenCoord; k++) {
            int i = icoordInner[k];
            int j = jcoordInner[k];
            double term1 = (s[(i+1)*Ny+j] - s[(i-1)*Ny+j])*(v[i*Ny+(j+1)] - v[i*Ny+(j-1)]);
            double term2 = (s[i*Ny+(j+1)] - s[i*Ny+(j-1)])*(v[(i+1)*Ny+j] - v[(i-1)*Ny+j]);
            double term3 = (v[i*Ny+(j+1)] - 2*v[i*Ny+j] + v[i*Ny+(j-1)])/dy/dy;
            double term4 = (v[(i+1)*Ny+j] - 2*v[i*Ny+j] + v[(i-1)*Ny+j])/dx/dx;
            temp[i*Ny+j] = dt/4/dx/dy * (term1 - term2) + dt/Re * (term3 + term4);
        }
        
        if (rank > 0) {
            // Send calculated local v matrix to root process
            MPI_Send(temp, Nxy, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else {
            // Receive each individual local segment and combine them on root for global
            for (int src = 1; src < nprocs; src++) {
                double* vTemp = new double[Nxy];
                MPI_Recv(vTemp, Nxy, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                cblas_daxpy(Nxy, 1.0, vTemp, 1, temp, 1);

		// Clean up
                delete[] vTemp;
            }
        }

    } else {
        // Serial code, calculate updated interior vorticity for entire matrix v
        for (unsigned int i = 1; i < (Nx-1); i++) {
            for (unsigned int j = 1; j < (Ny-1); j++) {
                double term1 = (s[(i+1)*Ny+j] - s[(i-1)*Ny+j])*(v[i*Ny+(j+1)] - v[i*Ny+(j-1)]);
                double term2 = (s[i*Ny+(j+1)] - s[i*Ny+(j-1)])*(v[(i+1)*Ny+j] - v[(i-1)*Ny+j]);
                double term3 = (v[i*Ny+(j+1)] - 2*v[i*Ny+j] + v[i*Ny+(j-1)])/dy/dy;
                double term4 = (v[(i+1)*Ny+j] - 2*v[i*Ny+j] + v[(i-1)*Ny+j])/dx/dx;
                temp[i*Ny+j] = dt/4/dx/dy * (term1 - term2) + dt/Re * (term3 + term4);
            }
        }
    }
    
    // Update current v at time t matrix with the new v matrix at time t + dt
    cblas_daxpy(Nxy, 1, temp, 1, v, 1);
    
    // Deallocate temp matrix
    delete[] temp;
}

/**
 * @brief Solution to the Poisson problem to compute the stream-function at time t + dt (equation 12)
 *        Integrate all the steps together to solve for stream-function and vorticity iteratively
 * 	  Create a new PoissonSolver instance and solve Poisson problem
 */
void LidDrivenCavity::Integrate() {
	
    // Calculate main diagonal and super diagonal values in the A matrix
    double seconddiag = 1 / (dx * dx);     		// diagonal values furthest away from the main diagonal
    double firstdiag = 1 / (dy * dy);	    		// diagonal values right above and below the main diagonal
    double diag = 2 * (firstdiag + seconddiag);         // main diagonal values
    
    // Poisson solver
    PoissonSolver* poissonSolver = new PoissonSolver();
    poissonSolver->SetVariables((int) Nx, (int) Ny, diag, firstdiag, seconddiag);
    
    // Parallel code
    if (nprocs > 1) {
        if (rank == 0) {
            // Generate Laplacian A matrix at the root process
            poissonSolver->GenerateScaLaPackAIntMatrix();
            ScaLaPackMatrix = poissonSolver->GetScaLaPackAIntMatrix();
            ScaLaPackMatrixNx = poissonSolver->GetScaLaPackAIntMatrixNx();
            ScaLaPackMatrixNy = poissonSolver->GetScaLaPackAIntMatrixNy();
        }
        // Send Laplacian A matrix to root process
        // Receive and set Laplacian A matrix from root process on other processes
        MPI_Bcast(&ScaLaPackMatrixNx, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ScaLaPackMatrixNy, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank > 0) {
            ScaLaPackMatrix = new double[ScaLaPackMatrixNx * ScaLaPackMatrixNy];
        }
        MPI_Bcast(ScaLaPackMatrix, ScaLaPackMatrixNx * ScaLaPackMatrixNy, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank > 0) {
            poissonSolver->SetScaLaPackAIntMatrix(ScaLaPackMatrix, ScaLaPackMatrixNx, ScaLaPackMatrixNy);
        }
        
        // Initialize BLACS for ScaLaPack
        poissonSolver->InitialiseScaLaPack(Px, Py);
        // Pre-factor A matrix for faster solve
        poissonSolver->PrefactorAIntMatrixParallel();
    } else {
        // Serial code, generate Laplacian A matrix
        poissonSolver->GenerateLaPackAIntMatrix();
        // Pre-factor A matrix for faster solve
        poissonSolver->PrefactorAIntMatrixSerial();
    }
    
    // Initialize current time variable
    double tCurrent = 0.0;
    do {
        SetInteriorVorticity();
        SetVorticityBCs();
        NextInteriorVorticity();
        
        // Broadcast assembled vorticity matrix to all other processes
        if (nprocs > 1) {
            MPI_Bcast(v, Nxy, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        
        // Set stream function and vorticity matrices for the Poisson solver
        poissonSolver->SetVectors(s, v);
        
        // Solve depending either serial or parallel
        if (nprocs > 1) {
            poissonSolver->SolveParallel();
        } else {
            poissonSolver->SolveSerial();
        }
        
        // Update current stream function at time t with next stream function at time t + dt
        if (rank == 0) {
            poissonSolver->Updatex(s);
        }
        
        // Broadcast updated stream function to all other processes
        if (nprocs > 1) {
            MPI_Bcast(s, Nxy, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        
        // Update current time to next time
        tCurrent += dt;
    } while (tCurrent < T);
}

/**
 * @brief Output stream function and vorticity matrices to txt file
 */
void LidDrivenCavity::GenerateData() {

    if (rank == 0) {
        // Create filename from parameters used
        string filename = "Plot/" +  to_string((int)Lx) + "_" + to_string((int)Ly) + "_" + to_string(Nx) + "_" + to_string(Ny) + "_" + to_string((int)Re) + "_data.txt";
        
        // Open file to write/overwrite
        ofstream outputfile(filename, ios::out | ios::trunc);
        if (outputfile.is_open()) {
            // Output the initial parameters
            outputfile << "Lx, Ly, Nx, Ny, dt, T, Re" << endl;
            outputfile << Lx << "," << Ly <<  "," << Nx <<  "," << Ny <<  "," << dt <<  "," << T <<  "," << Re << endl;
            
            // Export vorticity
            outputfile << endl << "vdata" << endl;
            for (unsigned int j = 0; j < Ny; j++) {
                for (unsigned int i = 0; i < Nx; i++) {
                    outputfile << v[j + i * Ny] << "    ";
                }
                outputfile << endl;
            }
            
            // Export stream function
            outputfile << endl << "sdata" << endl;
            for (unsigned int j = 0; j < Ny; j++) {
                for (unsigned int i = 0; i < Nx; i++) {
                    outputfile << s[j + i * Ny] << "    ";
                }
                outputfile << endl;
            }
        }
        cout << endl << endl;
    }
}
