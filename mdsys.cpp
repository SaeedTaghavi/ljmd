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
    std::string restart_file;

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
    n_atoms = stoi(options[0]);
    n_steps = stoi(options[9]);
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
    rx.reserve(n_atoms);
    ry.reserve(n_atoms);
    rz.reserve(n_atoms);
    vx.reserve(n_atoms);
    vy.reserve(n_atoms);
    vz.reserve(n_atoms);
    fx.reserve(n_atoms);
    fy.reserve(n_atoms);
    fz.reserve(n_atoms);

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
        if (n_line < n_atoms)
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

void Mdsys::calculate_kin_en()
{
    kin_en = 0.0;

    for (int i = 0; i < n_atoms; i++)
    {
        kin_en += 0.5 * MVSQ2E * mass * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    temp = 2.0 * kin_en / (3.0 * n_atoms - 3.0) / BOLTZ;    
}

void Mdsys::calculate_forces()
{
    double r,ffac;
    double r_x,r_y,r_z;

    // Zeros the potential energy and forces
    pot_en = 0.0;
    fill(fx.begin(), fx.end(), 0.0);
    fill(fy.begin(), fy.end(), 0.0);
    fill(fz.begin(), fz.end(), 0.0);

    for (int i = 0; i < n_atoms; i++)
    {
        for (int j = 0; j < n_atoms; j++)
        {
            if (i != j)
            {
                // Get distance between particle i and j
                r_x = pbc_helper(rx[i] - rx[j], 0.5 * box);
                r_y = pbc_helper(ry[i] - ry[j], 0.5 * box);
                r_z = pbc_helper(rz[i] - rz[j], 0.5 * box);
                r = std::sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

                // Compute force and energy if within cutoff
                if (r < rcut)
                {
                    ffac = -4.0 * epsilon * (-12.0 * std::pow(sigma / r, 12.0) / r + 6 * std::pow(sigma / r, 6.0) / r);  
                    pot_en += 0.5 * 4.0 * epsilon * (std::pow(sigma / r, 12.0) - std::pow(sigma / r, 6.0));

                    fx[i] += r_x / r * ffac;
                    fy[i] += r_y / r * ffac;
                    fz[i] += r_z / r * ffac;
                }   
            }   
        }   
    }
}

void Mdsys::calculate_vel_verlet()
{
    // First part: propagate velocities by half and positions by full step
    for (int i = 0; i < n_atoms; i++)
    {
        vx[i] += 0.5 * dt / MVSQ2E * fx[i] / mass;
        vy[i] += 0.5 * dt / MVSQ2E * fy[i] / mass;
        vz[i] += 0.5 * dt / MVSQ2E * fz[i] / mass;
        rx[i] += dt * vx[i];
        ry[i] += dt * vy[i];
        rz[i] += dt * vz[i];
    }

    this->calculate_forces();

    // Second part: propagate velocities by another half step
    for (int i = 0; i < n_atoms; i++)
    {
        vx[i] += 0.5 * dt / MVSQ2E * fx[i] / mass;
        vy[i] += 0.5 * dt / MVSQ2E * fy[i] / mass;
        vz[i] += 0.5 * dt / MVSQ2E * fz[i] / mass;
    }   
}

void Mdsys::output_helper(std::ofstream& energy_file, std::ofstream& trajectory_file)
{
    std::cout << "     " << nfi << "            " << temp << "            " << kin_en << "            " << pot_en << "            " << kin_en + pot_en << '\n';
    energy_file << "     " << nfi << "            " << temp << "            " << kin_en << "            " << pot_en << "            " << kin_en + pot_en << '\n';
    trajectory_file << n_atoms << '\n' << "nfi = " << nfi << " etot = " << kin_en + pot_en << '\n';
    for (int i = 0; i < n_atoms; i++)
    {
        trajectory_file << "Ar " << rx[i] << "     " << ry[i] << "     " << rz[i] << '\n';
    } 
}

double pbc_helper(double x, const double boxby2)
{
    while (x >  boxby2) x -= boxby2 + boxby2;
    while (x < -boxby2) x += boxby2 + boxby2;
    return x;
}
