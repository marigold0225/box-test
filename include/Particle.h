#ifndef PARTICLE_H
#define PARTICLE_H

#include <tuple> // For the velocity() function return type
#include <random>
#include <cmath>
class Simulation;

enum ParticleType
{
    PROTON,
    PION,
    DELTA
};

class Particle
{
public:
    static const double PROTON_MASS;
    static const double PION_MASS;
    static const double DELTA_MASS0;

    ParticleType type;
    int birthStep;
    double x, y, z;
    double px, py, pz;
    double energy;
    double mass;

    Particle(ParticleType type);

    // Move the particle based on its momentum and the time step
    void move(double dt);

    // Lorentz boost
    void lorentzBoost(double beta_x, double beta_y, double beta_z);

    // Set momentum and update energy
    void setMomentum(double new_px, double new_py, double new_pz);

    // Get the particle's velocity as a tuple
    std::tuple<double, double, double> velocity() const;

    // Set particle's momentum in a random direction with given magnitude
    void setRandomMomentumDirection(double magnitude, Simulation &sim);

    // Utility function to get mass based on particle type
    static double getMass(ParticleType type);
};

#endif
