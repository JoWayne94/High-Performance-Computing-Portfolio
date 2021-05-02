# Lid-Driven Cavity Project

This project implements a parallel numerical code to solve the vorticity-streamfunction formulation of the incompressible Naiver-Stokes equation in 2D using the finite difference method. In other words, it presents a solver for the lid-driven cavity problem.

# Pre-requisites

Please make sure the libraries below are installed on the system before running the code:

* BLAS - linear algebra package to perform basic vector and matrix multiplication
* BLACS - linear algebra oriented message passing interface (MPI)
* Boost - utilise the **program_options** library to obtain program options from command line
* LAPACK - linear algebra package to solve systems of simultaneous linear equations in serial
* MPI - message passing interface for parallel programming
* ScaLAPACK - linear algebra package to solve systems of simultaneous linear equations in parallel

# Compile the program

Make sure you are in the correct directory containing the **Makefile** and run

```shell
make
```

The **Makefile** will compile the program and produce the executable **Solve** file.

# Run test cases

One serial and one parallel test cases are declared in the **Makefile** for validation purposes. Their parameters are:

* Lx = Ly = 1.0
* Nx = Ny = 5
* dt = 0.0001
* T = 1.0
* Re = 100

# Serial test case

Run with single process (--np 1) and Px = Py = 1:

```shell
make testSerial
```

# Parallel test case

Run with 4 processes (--np 4) with Px = Py = 2, Px * Py = 4:

```shell
make testParallel
```

Output of the tests will be in the **Plot/** folder.

# Running for report output

An additional script named **./RunSolver** was implemented to run the solver with the given parameters required for the HPC Coursework report plots. To run the script, run the following command

```shell
./RunSolver.sh
```

This script will first compile the solver via the **`make`** command, then run the solver for the following parameters

* Lx = 1.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 1.0, Re = 100
* Lx = 1.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 1.0, Re = 400
* Lx = 1.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 5.0, Re = 1000
* Lx = 1.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 10.0, Re = 3200
* Lx = 1.0, Ly = 2.0, Nx = 161, Ny = 161, dt = 0.0005, T = 1.0, Re = 100
* Lx = 2.0, Ly = 1.0, Nx = 161, Ny = 161, dt = 0.0005, T = 1.0, Re = 100

All of these are run with 9 processes (--np 9) and Px = Py = 3. The data is saved in the **Data/** folder in the format `<Lx>_<Ly>_<Nx>_<Ny>_<Re>_data.txt`. Following this, the MATLAB script **GeneratePlots.m** is run, which generates the required plots and table data for the report. These plots are saved in the **Images/** folder. 



## Cleaning the folder

To clean the folder, i.e. remove the generated files, the following commands can be run.

### Remove generated object file (*.o)

```shell
make clean
```

### Remove generated plots (*.png)

```shell
make cleanImages
```

### Remove generated data file (*.txt)

```shell
make cleanData
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
--help				prints out the help message
--Lx <double>	specifies length of domain in the x-direction
--Ly <double>	specifies length of domain in the y-direction
--Nx <int>		specifies number of grid points in the x-direction
--Ny <int>		specifies number of grid points in the y-direction
--Px <int>		specifies number of partitions in the x-direction
--Py <int>		specifies number of partitions in the y-direction
--dt <double>	specifies the timestep
--T <double>	specifies the final time
--Re <double>	specifies the Reynolds number
```

Similarly, the data generated from running this will be saved in the **Data/** folder in the format `<Lx>_<Ly>_<Nx>_<Ny>_<Re>_data.txt`.



## Data verification

A seperate MATLAB implementation of the lid-driven cavity solver is included in the distribution. This file is named **LidDrivenCavitySolver.m** and can be used to test the output against that produced by the C++ implementation.



## Directory Structure

This explains the purpose of each file and folder in this folder

* [assignment.pdf](./assignment.pdf) - HPC coursework assignment brief
* [README.md](./README.md) - README markdown file
* [README.pdf](./README.pdf) - PDF version of the README file
* [repository.log](./repository.log) - Git repository log file
* [LidDrivenCavitySolver.cpp](./LidDrivenCavitySolver.cpp) - C++ file that accepts options for the solver and then sets up the solver for the lid-driven cavity problem.
* [LidDrivenCavity.cpp](./LidDrivenCavity.cpp) - C++ file that implements the serial and parallel LidDrivenCavity class member functions
* [Poisson2DSolver.cpp](./Poisson2DSolver.cpp) - C++ file that implements the serial and parallel Poisson2DSolver class member functions
* [LidDrivenCavity.h](./LidDrivenCavity.h) - Header file for the LidDrivenCavity class
* [Poisson2DSolver.h](./Poisson2DSolver.h) - Header file for the Poisson2DSolver class
* [GeneratePlots.m](./GeneratePlots.m) - MATLAB script to generate plots from the output files of the program
* [GenerateScalePlot.m](./GeneratePlot.m) - MATLAB script to generate scale plot from the output file (executionTime.txt) of ./TestSpeed.sh
* [LidDrivenCavitySolver.m](./LidDrivenCavitySolver.m) - MATLAB implementation of the serial lid driven cavity solver
* [Makefile](./Makefile) - File to help in the compilation of the executable, cleaning and testing
* [RunSolver.sh](./RunSolver.sh) - Shell script to run solver to produce required data for the report
* [TestSpeed.sh](./TestSpeed.sh) - Shell script to run solver for different number of processes and time it
* [Data/](./Data/) - Folder that contains the output data of the solver
* [Images/](./Images/) - Folder that contains the plots produced by MATLAB script
* [ReportFiles/](./ReportFiles/) - LaTEX file to compile report

The following are the output files from running different scripts

* **LidDrivenCavitySolver.o** - Object file of LidDrivenCavitySolver.cpp from running make
* **LidDrivenCavity.o** - Object file of LidDrivenCavity.cpp from running make
* **Poisson2DSolver.o** - Object file of Poisson2DSolver.cpp from running make
* **./Solve** - Executable linking LidDrivenCavitySolver.o, LidDrivenCavity.o and Poisson2DSolver.o
* **executionTime.txt** - Text file containing run times for different number of processes. Produced from running ./TestSpeed.sh
* **./Data/*.txt** - Text files containing the data produced from running ./Solve.
* **./Images/*.png** - Image files used in the report produced from running the MATLAB script GeneratePlots.m.
* **./ReportFiles/*** - Files produced when compiling Report.tex.



All files in this repository are required to run the code to generate data points that map velocity, vorticity and streamfunction values at different grid locations.   
Documentation is written within the source code as comments, describing the functions and purpose of each file.

**Input**: To vary initial and runtime conditions of the flow, please change the variable values within the ProgramInputs.cpp file.   
**Output**: If the code ran successfully, an Output.txt file containing velocity, vorticity, and streamfunction of all defined timesteps throughout the time interval specified OR once convergence is reached will appear in your directory.
______
Linux

Navigate to the directory containing all source files, then enter the following commands accordingly:
1.  cmake ./
2.  make
3.  mpirun -np *insert number of partitions in total* MainSolver
