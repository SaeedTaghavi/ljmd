/* 
 * Simple Lennard-Jones Potential Molecular Dynamics code with velocity verlet.
 * Units: 
 * Length = Angstrom; Mass = amu; Energy = kcal
 *
 * Molecular Dynamics System header file
 */

#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <charconv>
#include <sstream>
#include <iostream>

using namespace std;

// Physical constants:
// Boltzman constant in kcal/mol/K
constexpr double boltz {0.0019872067};

// m*v^2 in kcal/mol
constexpr double mvsq2e {2390.05736153349};

// Structure that holds information about the system
struct Mdsys 
{
    int natoms, nfi, nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    vector<double> rx, ry, rz, vx, vy, vz, fx, fy, fz;
};

// Structure initialization
void mdsys_init(Mdsys& sys, string filename);

// Calculates the kinectic energy of the system
// void kin_en(Mdsys& sys);

// // Calculates the forces on the system
// void forces(Mdsys& sys);

// // Calculates the velocity verlet
// void vel_verlet(Mdsys& sys);