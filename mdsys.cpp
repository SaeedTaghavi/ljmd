/* 
 * Simple Lennard-Jones Potential Molecular Dynamics code with velocity verlet.
 * Units: 
 * Length = Angstrom; Mass = amu; Energy = kcal
 *
 * Molecular Dynamics System file
 */

#include "mdsys.hpp"

Mdsys::Mdsys(std::string filename)
{
    // Supporting filenames
    std::string restart_file, trajectory_file, energy_file;

    // Output print frequency
    int print_freq;

    // Read parameters from file
    std::ifstream file {filename};

    if (!file)
    {
        throw std::runtime_error("Error opening the initial file! \n");
    }

    std::string line;
    std::string delimiter = " ";
    std::string token;

    // Fills vector with options
    std::array<std::string, 12> options;
    int index {0};
    while (getline(file, line) && index < 12) 
    {
        token = line.substr(0, line.find(delimiter));
        options[index] = token;
        ++index;
    }

    // Configure atributes with the parameters
    natoms = stoi(options[0]);
    nsteps = stoi(options[9]);
    mass = stod(options[1]);
    epsilon = stod(options[2]);
    sigma = stod(options[3]);
    rcut = stod(options[4]);
    box = stod(options[5]);
    dt = stod(options[10]); 

    // Supporting files
    restart_file = "data/" + options[6];
    trajectory_file = options[7];
    energy_file = options[8];

    // Print frequency
    print_freq = stoi(options[11]);

    // Reserve memory for all vectors
    rx.reserve(natoms);
    ry.reserve(natoms);
    rz.reserve(natoms);
    vx.reserve(natoms);
    vy.reserve(natoms);
    vz.reserve(natoms);
    fx.reserve(natoms);
    fy.reserve(natoms);
    fz.reserve(natoms);

    // Initialize force vector elements as zero
    fill(fx.begin(), fx.end(), 0.0);
    fill(fy.begin(), fy.end(), 0.0);
    fill(fz.begin(), fz.end(), 0.0);

    // Initialize nfi as zero
    nfi = 0;

    // Read position and velocity from file
    std::ifstream res_file {restart_file};

    if (!res_file)
    {
        throw std::runtime_error("Error opening the restart file! \n");
    }

    int n_line {0};
    double col1, col2, col3;
    while (getline(res_file, line)) 
    {
        std::istringstream ss(line);
        if (n_line < natoms)
        {
            ss >> col1 >> col2 >> col3;
            rx.push_back(col1);
            ry.push_back(col2);
            rz.push_back(col3);
        }
        else
        {
            ss >> col1 >> col2 >> col3;
            vx.push_back(col1);
            vy.push_back(col2);
            vz.push_back(col3);
        }
        ++n_line;
    }
}

Mdsys::~Mdsys()
{
    std::cout << "Deleted the system!" << '\n';
}
