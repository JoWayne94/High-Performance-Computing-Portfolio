# Lid-Driven Cavity Project

This project implements a parallel numerical code to solve the vorticity-streamfunction formulation of the incompressible Naiver-Stokes equation in 2D using the finite difference method. In other words, it presents a solver for the lid-driven cavity problem. Output of solver is used to generate data points that map velocity, vorticity and streamfunction values at different grid locations. Full report can be found [here](https://www.researchgate.net/publication/344228346_Lid-Driven_Cavity_Report).

## Pre-requisites

Please make sure the libraries below are installed on the system before running the code:

* BLAS - linear algebra package to perform basic vector and matrix multiplication
* BLACS - linear algebra oriented message passing interface (MPI)
* Boost - utilise the **program_options** library to obtain program options from command line
* LAPACK - linear algebra package to solve systems of simultaneous linear equations in serial
* MPI - message passing interface for parallel programming
* ScaLAPACK - linear algebra package to solve systems of simultaneous linear equations in parallel

____________________________________________
UNIX/Linux

## Compile the program

Navigate to the directory containing all source codes, then enter the following commands accordingly:

```shell
make
```

The **Makefile** will compile the program and produce the executable **Solve** file.

## Run test cases

One serial and one parallel test cases are declared in the **Makefile** for validation purposes. Their parameters are:

* Lx = Ly = 1.0
* Nx = Ny = 5
* dt = 0.0001
* T = 1.0
* Re = 100

### Serial test case

Run with single process (--np 1) and Px = Py = 1:

```shell
make testSerial
```

### Parallel test case

Run with 4 processes (--np 4) with Px = Py = 2, Px * Py = 4:

```shell
make testParallel
```

Output of the tests will be in the **Plot/** folder.

## Run Cases

Shell script named **RunSolver.sh** batch runs the solver with given parameters:

```shell
./RunSolver.sh
```

* Lx = 1.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 1.0, Re = 100
* Lx = 1.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 1.0, Re = 400
* Lx = 1.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 5.0, Re = 1000
* Lx = 1.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 10.0, Re = 3200
* Lx = 1.0, Ly = 2.0, Nx = 161, Ny = 161, dt = 0.0005, T = 1.0, Re = 100
* Lx = 2.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 1.0, Re = 100

16 processes (--np 16), Px = Py = 4 are used for the runs. Output is saved in the **Plot/** folder in a `<Lx>_<Ly>_<Nx>_<Ny>_<Re>_data.txt` format.

## Clean folder

To remove generated files, run the following commands:

### Remove generated object files (*.o)

```shell
make clean
```

### Remove generated data file (*.txt)

```shell
make cleanPlot
```

### Remove all generated files

```shell
make cleanAll
```

## Usage

Once the solver has been compiled, it can be used via the format

```shell
mpiexec --np <no. of procs> ./Solve [options] [args]
```

```
options:
--help         prints help message
--Lx <double>  length of domain in the x-direction
--Ly <double>  length of domain in the y-direction
--Nx <int>     number of grid points in the x-direction
--Ny <int>     number of grid points in the y-direction
--Px <int>     number of partitions in the x-direction
--Py <int>     number of partitions in the y-direction
--dt <double>  time step size
--T <double>   final time
--Re <double>  Reynolds number
```

## Directory Files

* [assignment.pdf](./assignment.pdf)                   - HPC coursework assignment brief
* [Lid-Driven-Cavity.ipynb](./Lid-Driven-Cavity.ipynb) - A seperate Python implementation of a serial lid-driven cavity solver
* [LidDrivenCavity.cpp](./LidDrivenCavity.cpp)         - Implement serial and parallel LidDrivenCavity class member functions
* [LidDrivenCavity.h](./LidDrivenCavity.h)             - Header file for the LidDrivenCavity class
* [LidDrivenCavitySolver.cpp](./LidDrivenCavitySolver.cpp) - Accept options for the solver and sets up the solver for the lid-driven cavity problem
* [Makefile](./Makefile)                   - Compilation of executable, cleaning, and test cases
* [./Plot/](./Plot/)                         - Folder containing output data of the solver to generate plots
* [PoissonSolver.cpp](./PoissonSolver.cpp) - Implement serial and parallel PoissonSolver class member functions
* [PoissonSolver.h](./PoissonSolver.h)     - Header file for the PoissonSolver class
* [README.md](./README.md)                 - README markdown file
* [repository.log](./repository.log)       - Git repository log file
* [RunSolver.sh](./RunSolver.sh)           - Shell script to run solver and produce required data
* [SpeedBatch.sh](./SpeedBatch.sh)         - Shell script to batch run solver for different number of processes and output time taken

Output files from running different scripts:

* **executionTimes.txt** - Text file containing run times for different number of processes. Output of ./SpeedBatch.sh
* **LidDrivenCavity.o** - Object file of LidDrivenCavity.cpp from running make
* **LidDrivenCavitySolver.o** - Object file of LidDrivenCavitySolver.cpp from running make
* **PoissonSolver.o** - Object file of PoissonSolver.cpp from running make
* **./Plot/*.txt** - Text files containing data produced from running ./Solve
* **./Solve** - Executable linking LidDrivenCavitySolver.o, LidDrivenCavity.o and PoissonSolver.o
