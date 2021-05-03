/**
 * @file LidDrivenCavity.h
 *
 * High-Performance Computing
 *
 * Header file of LidDrivenCavity.cpp
 *
 * Define all public classes and private variables
 *
 */
#ifndef LIDDRIVENCAVITY_H
#define LIDDRIVENCAVITY_H
#pragma once
#include <string>

using namespace std;

class LidDrivenCavity {

public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(unsigned int nx, unsigned int ny);
    void SetPartitionSize(unsigned int px, unsigned int py);
    void SetTimeStepSize(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetGridSpacing(double deltax, double deltay);
    void SetMPIConfig();
    
    // Main solver class functions
    void Initialise();
    void Integrate();
    void GenerateData();

private:
    // Lid velocity is 1
    const double U = 1.0;

    // Define parameters
    double dt;
    double T;
    unsigned int Nx;
    unsigned int Ny;
    unsigned int Px;
    unsigned int Py;
    double Lx;
    double Ly;
    double Re;
    double dx;
    double dy;
    double* v;
    double* s;
    
    unsigned int Nxy;
    unsigned int innerNx;
    unsigned int innerNy;
    unsigned int innerNxy;
    
    // Variables for MPI
    int rank;
    int nprocs;
    
    // Variables for individual partitions to work
    int arrlenCoord;
    int* icoordInner;
    int* jcoordInner;
    
    // Variables for poisson solver
    double* ScaLaPackMatrix = nullptr;
    int ScaLaPackMatrixNx;
    int ScaLaPackMatrixNy;
    
    // Private functions to solve vorticity and stream function
    void SetVorticityBCs();
    void SetInteriorVorticity();
    void NextInteriorVorticity();
};

#endif
