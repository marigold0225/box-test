#include "Particle.h"
#include "Box.h"
#include "Simulation.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>

std::map<std::string, double> readParametersFromFile(const std::string &filename)
{
    std::ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        std::cerr << "Error: Couldn't open the input file!" << std::endl;
        exit(1);
    }

    std::map<std::string, double> parameters;
    std::string line, key;
    double value;

    while (std::getline(inputFile, line))
    {
        std::istringstream iss(line);
        if (std::getline(iss, key, '=') && (iss >> value))
        {
            parameters[key] = value;
        }
    }
    inputFile.close();

    return parameters;
}

int main()
{
    std::map<std::string, double> parameters = readParametersFromFile("input.txt");

    // Create a simulation instance using the read parameters
    SimulationParameters config{
        static_cast<int>(parameters["numEvents"]),
        parameters["timestep"],
        parameters["totalTime"],
        static_cast<unsigned int>(parameters["seed"]),
        static_cast<int>(parameters["numProtons"]),
        static_cast<int>(parameters["numPions"]),
        parameters["kT"],
        parameters["boxSize"],
        parameters["desiredCellSize"],
        static_cast<int>(parameters["cdfResolution"]),
        parameters["maxP"]};

    Simulation simulation(config);

    // Print the initial state of the particles
    std::ofstream initFile("initial_state.csv");
    if (initFile.is_open())
    {
        initFile << "Type, Px, Py, Pz, X, Y, Z\n"; // header

        const auto &particles = simulation.getParticles();
        for (const Particle &particle : particles)
        {
            initFile << particle.type << ", "
                     << particle.px << ", "
                     << particle.py << ", "
                     << particle.pz << ", "
                     << particle.x << ", "
                     << particle.y << ", "
                     << particle.z << "\n";
        }
        initFile.close();
    }
    else
    {
        std::cerr << "Error: Unable to open file for writing initial state." << std::endl;
        return 1;
    }

    simulation.run();

    std::cout << "Simulation completed. Check output files for results." << std::endl;

    return 0;
}
