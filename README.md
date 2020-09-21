All files in this repository are required to run the code to generate data points that map velocity, vorticity and streamfunction values at different grid locations.
Documentation is written within the source code as comments, describing the functions and purpose of each file.

Input: To vary initial and runtime conditions of the flow, please change the variable values within the ProgramInputs.cpp file
Output: If the code ran successfully, an Output.txt file containing velocity, vorticity, and streamfunction of all defined timesteps throughout the time interval specified
        OR once convergence is reached, will appear in your directory
______
Linux

Navigate to the directory containing all source files, then enter the following commands accordingly:
1.  cmake ./
2.  make
3.  mpirun -np *insert number of partitions in total* MainSolver
