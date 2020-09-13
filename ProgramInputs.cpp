/**
 * @file ProgramInputs.cpp
 *
 * High-performance Computing
 *
 * Controls the input from command line, storing it in predefined variables
 */
#include <boost/program_options.hpp>
#include<cstdlib>
#include<iostream>

using namespace std;

namespace po = boost::program_options;

/**
 * @brief Read variables from vm or command line
 */
void Inputs(po::variables_map &vm, double &dt, double &T, int &Nx, int &Ny, double &Lx, double &Ly, double &Re, int &Px, int &Py) {
    dt = vm["dt"].as<double>();
    T =  vm["T"] .as<double>();
    Nx = vm["Nx"].as<int>();
    Ny = vm["Ny"].as<int>();
    Lx = vm["Lx"].as<double>();
    Ly = vm["Ly"].as<double>();
    Re = vm["Re"].as<double>();
    Px = vm["Px"].as<int>();
    Py = vm["Py"].as<int>();
}

/**
 * @brief Add input options with appropriate default values to run the code.
 */
bool InputStatus(int argc, char* argv[], po::variables_map &vm) {
   try
   {
    po::options_description desc("Default input options.");

    desc.add_options()
        ("help", "Help message.")
        ("Lx", po::value<double>()-> default_value(1.0), "Length of domain in x-direction.")
        ("Ly", po::value<double>()-> default_value(1.0), "Length of domain in y-direction.")
        ("Nx", po::value<int>()   -> default_value(161), "Number of grid points in x-direction.")
        ("Ny", po::value<int>()   -> default_value(161), "Number of grid points in y-direction.")
        ("Px", po::value<int>()   -> default_value(4),   "Number of partitions in x-direction (parallel)")
        ("Py", po::value<int>()   -> default_value(4),   "Number of partitions in y-direction (parallel)")
        ("dt", po::value<double>()-> default_value(0.00976), "Time step size.")
        ("T",  po::value<double>()-> default_value(999.0), "Final time.")
        ("Re", po::value<double>()-> default_value(1000.0), "Reynolds number.");

    // Parse command-line arguments and store in buffer vm
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
   }
   catch(exception const &e)
   {
    cout << e.what() << endl;
    return false;
   }
   return true;
}