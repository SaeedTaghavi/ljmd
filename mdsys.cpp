/* 
 * Simple Lennard-Jones Potential Molecular Dynamics code with velocity verlet.
 * Units: 
 * Length = Angstrom; Mass = amu; Energy = kcal
 *
 * Molecular Dynamics System file
 */

#include "mdsys.hpp"

// Structure initialization
void mdsys_init(Mdsys& sys, string filename) 
{
    // Supporting filenames
    string restart_file, trajectory_file, energy_file;

    // Output print frequency
    int print_freq;

    // Read parameters from file
    ifstream file {filename};

    if (!file)
    {
        throw runtime_error("Error opening the initial file! \n");
    }
    
    string line;
    string delimiter = " ";
    string token;

    array<string, 12> options;
    int index {0};
    while (getline(file, line) && index < 12) 
    {
        token = line.substr(0, line.find(delimiter));
        options[index] = token;
        ++index;
    }

    // Configure atributes with the parameters
    from_chars(options[0].data(), options[0].data() + options[0].size(), sys.natoms);
    from_chars(options[9].data(), options[9].data() + options[9].size(), sys.nsteps);
    sys.mass = stod(options[1]);
    sys.epsilon = stod(options[2]);
    sys.sigma = stod(options[3]);
    sys.rcut = stod(options[4]);
    sys.box = stod(options[5]);
    sys.dt = stod(options[10]); 

    // Supporting files
    restart_file = options[6];
    trajectory_file = options[7];
    energy_file = options[8];

    // Print frequency
    from_chars(options[11].data(), options[11].data() + options[11].size(), print_freq);

    // Reserve memory for all vectors
    sys.rx.reserve(sys.natoms);
    sys.ry.reserve(sys.natoms);
    sys.rz.reserve(sys.natoms);
    sys.vx.reserve(sys.natoms);
    sys.vy.reserve(sys.natoms);
    sys.vz.reserve(sys.natoms);
    sys.fx.reserve(sys.natoms);
    sys.fy.reserve(sys.natoms);
    sys.fz.reserve(sys.natoms);

    // Initialize force vector elements as zero
    fill(sys.rx.begin(), sys.fx.end(), 0.0);
    fill(sys.ry.begin(), sys.fx.end(), 0.0);
    fill(sys.rz.begin(), sys.fx.end(), 0.0);
    fill(sys.vx.begin(), sys.fx.end(), 0.0);
    fill(sys.vy.begin(), sys.fx.end(), 0.0);
    fill(sys.vz.begin(), sys.fx.end(), 0.0);
    fill(sys.fx.begin(), sys.fx.end(), 0.0);
    fill(sys.fy.begin(), sys.fx.end(), 0.0);
    fill(sys.fz.begin(), sys.fx.end(), 0.0);

    // Read position and velocity from file
    ifstream res_file {restart_file};

    if (!res_file)
    {
        throw runtime_error("Error opening the restart file! \n");
    }

    int n_line {0};
    double col1, col2, col3;
    cout << n_line << '\n';
    while (getline(res_file, line)) 
    {
        istringstream ss(line);
        cout << line << '\n';
        if (n_line < sys.natoms)
        {
            ss >> col1 >> col2 >> col3;
            cout << n_line << '\n';
            sys.rx.at(n_line) = col1;
            sys.ry.at(n_line) = col2;
            sys.rz.at(n_line) = col3;
        }
        else
        {
            ss >> col1 >> col2 >> col3;
            cout << n_line << '\n';
            sys.vx.at(n_line - sys.natoms) = col1;
            sys.vy.at(n_line - sys.natoms) = col2;
            sys.vz.at(n_line - sys.natoms) = col3;
        }
        ++n_line;
    }

    
}
