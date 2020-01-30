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

// Physical constants:
// Boltzman constant in kcal/mol/K
const double BOLTZ {0.0019872067};

// m*v^2 in kcal/mol
const double MVSQ2E {2390.05736153349};

// Class that holds information about the system
class Mdsys
{
private:
    int natoms, nfi, nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    std::vector<double> rx, ry, rz, vx, vy, vz, fx, fy, fz;
public:
    Mdsys(std::string filename);

    ~Mdsys();

    // Calculates the kinectic energy of the system
    //void kin_en();

    // Calculates the forces on the system
    // void forces();

    // Calculates the velocity verlet
    // void vel_verlet();
};
