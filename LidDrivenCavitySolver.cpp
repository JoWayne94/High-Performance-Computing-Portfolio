/**
 * @file LidDrivenCavitySolver.cpp
 *
 * High-Performance Computing
 *
 * Solution to the Coursework Assignment
 *
 * Main file of the solution
 * Controls the input from command line, storing it in pre-defined variables
 * Accept options for the solver and sets up the solver for the Lid-Driven Cavity problem
 */
#include "LidDrivenCavity.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <exception>
#include <iomanip>
#include <cmath>
#include <mpi.h>
#include <cstdlib>
#include <iterator>
using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv) {

    // Set up command line options and defaults. Add input options with appropriate default values to run the code
    po::options_description opts("Solves the lid-driven cavity problem using given algorithms");
	opts.add_options()
	("help", "Prints help messege")
	("Lx", po::value<double>()      ->default_value(1.0), "Length of domain in x-direction. Default = 1.0")
        ("Ly", po::value<double>()      ->default_value(1.0), "Length of domain in y-direction. Default = 1.0")
        ("Nx", po::value<unsigned int>()->default_value(161), "Number of grid points in x-direction. Default = 161")
        ("Ny", po::value<unsigned int>()->default_value(161), "Number of grid points in y-direction. Default = 161")
        ("Px", po::value<unsigned int>()->default_value(1),   "Number of partitions in x-direction. Default = 1 (serial)")
        ("Py", po::value<unsigned int>()->default_value(1),   "Number of partitions in y-direction. Default = 1 (serial)")
        ("dt", po::value<double>()      ->default_value(0.0001), "Time step size. Default = 1E-4")
        ("T",  po::value<double>()      ->default_value(1),      "Final time. Default = 1.0 second(s)")
        ("Re", po::value<double>()      ->default_value(100),    "Reynolds number. Default = 100");
    
    // Tell boost to parse command-line arguments using list of possible options and generate a map of options and values
    po::variables_map vm;
    // Parse command-line arguments and store in buffer vm
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);
    
    // Check if user used the "--help" option and print usage
    if (vm.count("help")) {
        cout << "Solves the lid-driven cavity problem for the given inputs" << endl;
        cout << opts << endl;
        return 0;
    }
    
    // Read variables from vm or command line. Extract options and save as variables
    const double Lx = vm["Lx"].as<double>();
    const double Ly = vm["Ly"].as<double>();
    const double Nx = vm["Nx"].as<unsigned int>();
    const double Ny = vm["Ny"].as<unsigned int>();
    const double Px = vm["Px"].as<unsigned int>();
    const double Py = vm["Py"].as<unsigned int>();
    const double dt = vm["dt"].as<double>();
    const double T  = vm["T"].as<double>();
    const double Re = vm["Re"].as<double>();
    
    // Check if partitions match MPI. Initialize variables for MPI to return
    int rank = 0;
    int size = 0;

    // Initialise MPI
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        cout << "Error: Failed to initialise MPI" << endl;
        return -1;
    }

    // Get comm rank and size of each process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check if number of processes match Nx * Ny
    if ((Px * Py) != (unsigned int) size) {
        if (rank == 0) {
            cout << "Error: Number of processes must equal to Nx * Ny" << endl;
        }
        MPI_Finalize();
        return 0;
    }
    
    // Display chosen parameter values
    if (rank == 0) {
        cout << "Selected length of domain in x-direction = " << Lx << endl;
        cout << "Selected length of domain in y-direction = " << Ly << endl;
        cout << "Selected number of grid points in x-direction = " << Nx << endl;
        cout << "Selected number of grid points in y-direction = " << Ny << endl;
        cout << "Selected number of partitions in x-direction = " << Px << endl;
        cout << "Selected number of partitions in y-direction = " << Py << endl;
        cout << "Selected time step size = " << dt << endl;
        cout << "Selected final time = " << T << endl;
        cout << "Selected Reynolds number = " << Re << endl;
        cout << endl;
    }
    
    // Perform checks on input variables. Check if grid length is valid (positive)
    if ((Lx <= 0) || (Ly <= 0)) {
        if (rank == 0) {
            cout << "Error: Minimum length of the domain must be > 0" << endl;
        }
        return -1;
    }
    
    // Check if discretisation is valid (more than 3 for second order central difference scheme)
    if (Nx < 3 || Ny < 3) {
        if (rank == 0) {
            cout << "Error: Minimum number of discretisation in x and y direction is 3. Please pick a value >= 3" << endl;
        }
        return -1;
    }
    
    // Check if time step size, final time, and Reynolds number are positive
    if ((T <= 0) || (dt <= 0) || (Re <= 0)) {
        if (rank == 0) {
            cout << "Error: Minimum values of dt, T, and/or Re must be > 0" << endl;
        }
        return -1;
    }
    
    // Caclulate grid spacing
    const double dx = Lx / ((double) Nx - 1);
    const double dy = Ly / ((double) Ny - 1);
    
    // Display grid spacing
    if (rank == 0) {
        cout << "Grid spacing in the x-direction is " << dx << endl;
        cout << "Grid spacing in the y-direction is " << dy << endl;
        cout << endl;
    }
    
    // Check if chosen time step meets the stability condition
    if (dt >= (Re * dx * dy / 4)) {
        if (rank == 0) {
            cout << "Re * dx * dy / 4 = " << Re * dx * dy / 4;
            cout << "Error: Time step size chosen is too large. Please choose a dt < Re*dt*dx/4 to meet stability requirements" << endl;
        }
        return -1;
    }
    
    // Creates new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
    
    // Configure solver
    solver->SetDomainSize(Lx, Ly);
    solver->SetGridSize(Nx, Ny);
    solver->SetPartitionSize(Px, Py);
    solver->SetTimeStepSize(dt);
    solver->SetFinalTime(T);
    solver->SetReynoldsNumber(Re);
    solver->SetGridSpacing(dx,dy);
    solver->SetMPIConfig();

    // Run solver
    solver->Initialise();
    solver->Integrate();
    solver->GenerateData();
    
    // Finalise MPI.
    MPI_Finalize();
    
    return 0;
}
