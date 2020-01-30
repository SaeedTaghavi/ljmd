/* 
 * Simple Lennard-Jones Potential Molecular Dynamics code with velocity verlet.
 * Units: 
 * Length = Angstrom; Mass = amu; Energy = kcal
 *
 */

#include <iostream>
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

    return 0;
}
