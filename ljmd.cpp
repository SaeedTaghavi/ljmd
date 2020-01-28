/* 
 * Simple Lennard-Jones Potential Molecular Dynamics code with velocity verlet.
 * Units: 
 * Length = Angstrom; Mass = amu; Energy = kcal
 *
 * Translation to C++ from the original in C by Axel Kohlmeyer
 */

#include <iostream>

#include "mdsys.hpp"

int main(int argc, char const *argv[])
{
    cout << "LJMD - C++ \n\n";

    if (argc != 2)
    {
        cout << "USAGE: ./ljmd [filename] \n";
        cout << "Aborting due to usage error \n";
        
        // EINVAL: Invalid argument
        return 22; 
    }
    
    Mdsys sys;
    mdsys_init(sys, argv[1]);

    cout << "Posição 5 do rx: \n";
    cout << sys.rx.at(4) << '\n';
    cout << "Posição 10 do fz: \n";
    cout << sys.fz.at(9) << '\n';
    cout << "Posição 2 do vy: \n";
    cout << sys.vy.at(1) << '\n';

    return 0;
}
