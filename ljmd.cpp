/* 
 * Simple Lennard-Jones Potential Molecular Dynamics code with velocity verlet.
 * Units: 
 * Length = Angstrom; Mass = amu; Energy = kcal
 *
 */

#include <memory>

#include "mdsys.hpp"

int main(int argc, char const *argv[])
{
    std::cout << "LJMD - C++ \n\n";

    if (argc != 2)
    {
        std::cout << "USAGE: ./ljmd [filename] \n";
        std::cout << "Aborting due to usage error \n";
        
        // EINVAL: Invalid argument
        return 22; 
    }
    
    // sys is an unique pointer to a Mdsys object
    std::unique_ptr<Mdsys> sys = std::make_unique<Mdsys>(argv[1]);

    // Initialize forces and energy
    sys->calculate_forces();
    sys->calculate_kin_en();

    // Open trajectory and energy files
    std::ofstream en_file {sys->energy_file};
    std::ofstream traj_file {sys->trajectory_file};

    std::cout << "Starting simulation with " << sys->n_atoms << " atoms for " << sys->n_steps << " steps.\n";
    std::cout << "     NFI            TEMP            EKIN                 EPOT              ETOT\n";

    sys->output_helper(en_file, traj_file);

    // Main MD Loop
    for(sys->nfi = 1; sys->nfi <= sys->n_steps; ++sys->nfi) {

        // Write output, if requested
        //if ((sys->nfi % sys->print_freq) == 0)
            sys->output_helper(en_file, traj_file);

        // Propagate system and recompute energies
        sys->calculate_vel_verlet();
        sys->calculate_kin_en();

    }

    std::cout << "Simulation completed!" << '\n';

    return 0;
}
