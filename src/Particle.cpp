#include "Particle.h"
#include "Simulation.h"

const double Particle::PROTON_MASS = 0.938;
const double Particle::PION_MASS = 0.135;
const double Particle::DELTA_MASS0 = 1.232;

Particle::Particle(ParticleType type) : type(type),
                                        x(0.0), y(0.0), z(0.0),
                                        px(0.0), py(0.0), pz(0.0),
                                        energy(0.0), birthStep(-1)
{
    mass = getMass(type);
}

double Particle::getMass(ParticleType type)
{
    switch (type)
    {
    case PROTON:
        return PROTON_MASS;
    case PION:
        return PION_MASS;
    case DELTA:
        return DELTA_MASS0;
    default:
        return 0.0;
    }
}

void Particle::move(double dt)
{
    double vx, vy, vz;
    std::tie(vx, vy, vz) = velocity();
    x += vx * dt;
    y += vy * dt;
    z += vz * dt;
}

void Particle::lorentzBoost(double beta_x, double beta_y, double beta_z)
{
    double beta2 = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z;
    double gamma = 1.0 / std::sqrt(1.0 - beta2);

    // Compute the dot product of beta and p for later use
    double beta_dot_p = beta_x * px + beta_y * py + beta_z * pz;

    double E_prime = gamma * (energy + beta_dot_p);
    double gamma2_over_one_plus_gamma = gamma * gamma / (1.0 + gamma);

    px += gamma * beta_x * energy - gamma2_over_one_plus_gamma * beta_dot_p * beta_x;
    py += gamma * beta_y * energy - gamma2_over_one_plus_gamma * beta_dot_p * beta_y;
    pz += gamma * beta_z * energy - gamma2_over_one_plus_gamma * beta_dot_p * beta_z;
    energy = E_prime;
}

void Particle::setMomentum(double new_px, double new_py, double new_pz)
{
    px = new_px;
    py = new_py;
    pz = new_pz;
    energy = std::sqrt(px * px + py * py + pz * pz + mass * mass);
}

std::tuple<double, double, double> Particle::velocity() const
{
    return std::make_tuple(px / energy, py / energy, pz / energy);
}

void Particle::setRandomMomentumDirection(double magnitude, Simulation &sim)
{
    double theta = 2 * acos(-1) * sim.getRandomNumber();
    double phi = std::acos(2 * sim.getRandomNumber() - 1);

    px = magnitude * std::sin(phi) * std::cos(theta);
    py = magnitude * std::sin(phi) * std::sin(theta);
    pz = magnitude * std::cos(phi);
    energy = std::sqrt(px * px + py * py + pz * pz + mass * mass);
}
