/**
 * @file ProgramInputs.h
 *
 * High-performance Computing
 *
 * Header file of ProgramInputs.cpp
 *
 * Define functions to read input variables
 */
#include <boost/program_options.hpp>

#pragma once
using namespace std;

namespace po = boost::program_options;

bool InputStatus(int argc, char* argv[], po::variables_map &vm);
void Inputs(po::variables_map &vm, double &dt, double &T, int &Nx, int &Ny, double &Lx, double &Ly, double &Re, int &Px, int &Py);