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
#include <cmath>
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
    double dt, mass, epsilon, sigma, box, rcut;
    double kin_en, pot_en, temp;
    std::vector<double> rx, ry, rz, vx, vy, vz, fx, fy, fz;

public:
    std::string trajectory_file, energy_file;
    int nfi, n_atoms, n_steps, print_freq;

public:
    Mdsys(std::string filename);

    ~Mdsys();

    // Calculates the kinectic energy of the system
    void calculate_kin_en();

    // Calculates the forces on the system
    void calculate_forces();

    // Calculates the velocity verlet
    void calculate_vel_verlet();

    // Prints and writes the results to files
    void output_helper(std::ofstream& energy_file, std::ofstream& trajectory_file);
};

double pbc_helper(double x, const double boxby2);
