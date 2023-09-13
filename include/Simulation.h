#ifndef SIMULATION_H
#define SIMULATION_H

#include "Particle.h"
#include "Box.h"
#include <vector>
#include <map>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <iostream>
#include <unordered_set>

struct SimulationParameters
{
    int numEvents;
    double timestep;
    double totalTime;
    unsigned int customSeed;
    int numProtons;
    int numPions;
    double kT;
    double boxSize;
    double desiredCellSize;
    int cdfResolution;
    double maxP;
};

class Simulation
{
public:
    explicit Simulation(const SimulationParameters &params);
    // random number generator
    void setSeed(unsigned int newSeed);
    unsigned int getSeed() const;
    double getRandomNumber();
    std::mt19937 &getRandomGenerator();

    // main
    void run();
    void processInteractions(int currentStep);
    // initialize function
    void initializeParticles();
    double gauss_integral(double down, double up, double mass);
    void computeCDF(double mass);
    double sampleBoltzmann(double mass);

    // collision part
    bool collisionOccurs(const Particle &particle1, const Particle &particle2, int currentStep);
    double collisionProbability(const Particle &particle1, const Particle &particle2);
    Particle createDeltaFromCollision(const Particle &p1, const Particle &p2, int currentStep);
    // ElasticCollision
    //    bool elasticCollisionOccurs(const Particle &particle1, const Particle &particle2);
    //    void handleElasticCollision(Particle &particle1, Particle &particle2);
    // decay part
    bool decayOccurs(const Particle &deltaParticle);
    double decayProbability(const Particle &particle);
    std::pair<Particle, Particle> createProtonAndPionFromDecay(const Particle &particle, int currentStep);

    // other
    double adjustForPeriodicBoundary(double coord1, double coord2) const;
    void updateParticlePositions();
    void outputEventHeader(int eventNumber);
    void outputCurrentState(int step);
    std::vector<Particle> &getParticles() { return box.getParticles(); }

private:
    int numEvents;
    double dt;
    double totalTime;
    int numProtons;
    int numPions;
    double kT;
    double boxSize;
    double desiredCellSize;
    int CDF_RESOLUTION;
    double MAX_P;
    Box box;
    unsigned int seed;

    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;
    std::vector<double> cdf_values;
    std::ofstream outFile;

    // Gauss integral constants
    const std::vector<double> weights = {0.3626837833783620, 0.3626837833783620,
                                         0.3137066458778873, 0.3137066458778873,
                                         0.2223810344533745, 0.2223810344533745,
                                         0.1012285362903763, 0.1012285362903763};
    const std::vector<double> nodes = {-0.1834346424956498, 0.1834346424956498,
                                       -0.5255324099163290, 0.5255324099163290,
                                       -0.7966664774136267, 0.7966664774136267,
                                       -0.9602898564975363, 0.9602898564975363};
};

#endif // SIMULATION_H
