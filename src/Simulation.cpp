#include "Simulation.h"
#include "Particle.h"
#include "Box.h"

// initialize Simulation function
Simulation::Simulation(const SimulationParameters &params)
    : numEvents(params.numEvents),
      dt(params.timestep),
      totalTime(params.totalTime),
      numProtons(params.numProtons),
      numPions(params.numPions),
      kT(params.kT),
      boxSize(params.boxSize),
      desiredCellSize(params.desiredCellSize),
      CDF_RESOLUTION(params.cdfResolution),
      MAX_P(params.maxP),
      dist(0.0, 1.0)
{
    seed = params.customSeed;
    gen.seed(seed);
    box.setSize(params.boxSize);
    box.initializeGrid(params.desiredCellSize);
    initializeParticles();
}
// seed and randomNumber
void Simulation::setSeed(unsigned int newSeed)
{
    seed = newSeed;
    gen.seed(seed);
}

unsigned int Simulation::getSeed() const
{
    return seed;
}

double Simulation::getRandomNumber()
{
    return dist(gen);
}

std::mt19937 &Simulation::getRandomGenerator()
{
    return gen;
}
// initialization of particles

double Simulation::gauss_integral(double down, double up, double mass)
{
    double result = 0.0;
    double t, f_val;

    for (int i = 0; i < nodes.size(); i++)
    {
        t = (up - down) * nodes[i] / 2 + (up + down) / 2;
        f_val = 4 * acos(-1) * t * t * std::exp(-std::sqrt(t * t + mass * mass) / kT) / (pow(2 * acos(-1), 3));
        result += weights[i] * f_val;
    }

    result *= (up - down) / 2;
    return result;
}

void Simulation::computeCDF(double mass)
{
    cdf_values.clear();
    double dp = MAX_P / CDF_RESOLUTION;
    double current_sum = 0.0;

    for (int i = 0; i < CDF_RESOLUTION; i++)
    {
        current_sum += gauss_integral(i * dp, (i + 1) * dp, mass);
        cdf_values.push_back(current_sum);
    }
}

double Simulation::sampleBoltzmann(double mass)
{
    double random_value = dist(gen) * cdf_values.back(); // Normalize by the last value

    // Use binary search to find the closest CDF value
    auto iter = std::lower_bound(cdf_values.begin(), cdf_values.end(), random_value);

    // Convert iterator to index
    int index = std::distance(cdf_values.begin(), iter);

    // Convert index to momentum value
    return (index + 0.5) * (MAX_P / CDF_RESOLUTION); // +0.5 to get the center of the bin
}

void Simulation::initializeParticles()
{
    box.clearParticles();
    // Compute CDF for protons
    computeCDF(Particle::getMass(PROTON));

    for (int i = 0; i < numProtons; ++i)
    {
        Particle proton(PROTON);
        proton.x = boxSize * dist(gen);
        proton.y = boxSize * dist(gen);
        proton.z = boxSize * dist(gen);

        double p = sampleBoltzmann(Particle::getMass(PROTON));
        proton.setRandomMomentumDirection(p, *this);
        proton.energy = sqrt(proton.px * proton.px + proton.py * proton.py +
                             proton.pz * proton.pz + Particle::getMass(PROTON) * Particle::getMass(PROTON));
        box.addParticle(proton);
    }

    // Compute CDF for pions - only if different from protons
    if (Particle::getMass(PION) != Particle::getMass(PROTON))
        computeCDF(Particle::getMass(PION));

    for (int i = 0; i < numPions; ++i)
    {
        Particle pion(PION);
        pion.x = boxSize * dist(gen);
        pion.y = boxSize * dist(gen);
        pion.z = boxSize * dist(gen);

        double p = sampleBoltzmann(Particle::getMass(PION));
        pion.setRandomMomentumDirection(p, *this);
        pion.energy = sqrt(pion.px * pion.px + pion.py * pion.py +
                           pion.pz * pion.pz + Particle::getMass(PION) * Particle::getMass(PION));
        box.addParticle(pion);
    }

    //    std::cout << "Initialization completed." << std::endl;
}

// collision process
double Simulation::collisionProbability(const Particle &particle1, const Particle &particle2)
{
    const double M_proton = Particle::getMass(PROTON);
    const double M_pion = Particle::getMass(PION);
    const double M_delta0 = Particle::getMass(DELTA);

    double totalEnergy = particle1.energy + particle2.energy;
    double totalPx = particle1.px + particle2.px;
    double totalPy = particle1.py + particle2.py;
    double totalPz = particle1.pz + particle2.pz;

    double scm = totalEnergy * totalEnergy - (totalPx * totalPx + totalPy * totalPy + totalPz * totalPz);
    double p_lab = (scm - M_proton * M_proton - M_pion * M_pion) / 2 / M_proton;
    p_lab = std::sqrt(pow(p_lab, 2) - M_pion * M_proton);
    double q = sqrt((scm - pow((M_proton + M_pion), 2)) * (scm - pow((M_proton - M_pion), 2))) / (2 * sqrt(scm));

    double Gamma = 0.47 * q * q * q / (M_pion * M_pion + 0.6 * q * q);
    double A = 4 * M_delta0 * M_delta0 * Gamma / (pow((scm - M_delta0 * M_delta0), 2) + M_delta0 * M_delta0 * Gamma * Gamma) / 0.948;

    auto [vx1, vy1, vz1] = particle1.velocity();
    auto [vx2, vy2, vz2] = particle2.velocity();

    double v_rel_x = (vx1 - vx2);
    double v_rel_y = (vy1 - vy2);
    double v_rel_z = (vz1 - vz2);
    double v = sqrt(v_rel_x * v_rel_x + v_rel_y * v_rel_y + v_rel_z * v_rel_z);

    double cross_sections = 2 * acos(-1) * Gamma * A / q / q;

    double dV = 27 * box.getCellVolume();

    double P = 0.197 * 0.197 * cross_sections * v * dt / dV;

    if (P > 1.0 || P < 0.0)
    {
        std::cerr << "Error: Invalid probability value for collision: " << P << std::endl;
        exit(1); // or handle the error in another way
    }

    return P;
}

bool Simulation::collisionOccurs(const Particle &particle1, const Particle &particle2, int currentStep)
{

    if (particle1.birthStep == currentStep && particle2.birthStep == currentStep)
    {
        //        std::cout << "Collision prevented due to same birthStep: " << particle1.birthStep << std::endl;
        return false;
    }

    double collisionProb = collisionProbability(particle1, particle2);

    return dist(gen) < collisionProb;
}

double Simulation::adjustForPeriodicBoundary(double coord1, double coord2) const
{

    double boxSize = box.getSize();
    double gridSize = box.getGridSize();
    double maxDist = 1.5 * gridSize; // Maximum distance between two neighboring particles considering 3x3x3 grid

    double diff = coord2 - coord1;

    if (std::abs(diff) <= maxDist)
    {
        // Particles are close enough, no need to adjust for periodic boundary
        return (coord1 + coord2) / 2.0;
    }
    else
    {
        // Here, we further check which particle is closer to the right boundary
        if (boxSize - std::max(coord1, coord2) < std::min(coord1, coord2))
        {
            // The particle closer to the right boundary gets adjusted for periodic boundary
            return (coord1 + coord2 - boxSize) / 2.0;
        }
        else
        {
            // The particle closer to the left boundary gets adjusted for periodic boundary
            return (coord1 + coord2 + boxSize) / 2.0;
        }
    }
}

Particle Simulation::createDeltaFromCollision(const Particle &p1, const Particle &p2, int currentStep)
{

    Particle delta(DELTA);

    // Calculate energy and momentum for the new DELTA particle
    double E = p1.energy + p2.energy;
    double Px = p1.px + p2.px;
    double Py = p1.py + p2.py;
    double Pz = p1.pz + p2.pz;

    delta.energy = E;
    delta.px = Px;
    delta.py = Py;
    delta.pz = Pz;
    delta.mass = sqrt(E * E - Px * Px - Py * Py - Pz * Pz);

    delta.x = adjustForPeriodicBoundary(p1.x, p2.x);
    delta.y = adjustForPeriodicBoundary(p1.y, p2.y);
    delta.z = adjustForPeriodicBoundary(p1.z, p2.z);
    delta.birthStep = currentStep;
    return delta;
}

// handle collision
// bool Simulation::elasticCollisionOccurs(const Particle &particle1, const Particle &particle2)
//{
//}
//
// void Simulation::handleElasticCollision(Particle &particle1, Particle &particle2)
//{
//}

// decay process --------------------------------
double Simulation::decayProbability(const Particle &particle)
{

    const double M_delta0 = Particle::getMass(DELTA);
    const double M_pion = Particle::getMass(PION);
    const double M_proton = Particle::getMass(PROTON);
    double scm = M_delta0 * M_delta0;

    double q = sqrt((scm - pow((M_proton + M_pion), 2)) * (scm - pow((M_proton - M_pion), 2))) / (2 * sqrt(scm));

    double Gamma = 0.47 * q * q * q / (M_pion * M_pion + 0.6 * q * q);

    double gamma = particle.energy / M_delta0;

    double P = Gamma * dt / gamma / 0.197;

    if (P > 1.0 || P < 0.0)
    {
        std::cerr << "Error: Invalid probability value for decay: " << P << std::endl;
        exit(1); // or handle the error in another way
    }

    return P;
}

bool Simulation::decayOccurs(const Particle &deltaParticle)
{
    // Ensure the particle is of type DELTA
    if (deltaParticle.type != DELTA)
        return false;

    double decayProb = decayProbability(deltaParticle);

    // Return true if the random value is less than the decay probability
    return dist(gen) < decayProb;
}

std::pair<Particle, Particle> Simulation::createProtonAndPionFromDecay(const Particle &deltaParticle, int currentStep)
{
    Particle proton(ParticleType::PROTON);
    Particle pion(ParticleType::PION);

    // 1. In the rest frame of the DELTA particle, the momentum of PROTON and PION are opposite.
    // Calculate the magnitude of momentum based on energy conservation.
    double deltaRestEnergy = deltaParticle.mass; // This is the energy of the DELTA in its rest frame
    double momentumMagnitude = sqrt((pow(deltaRestEnergy, 2) - pow(Particle::getMass(PROTON) + Particle::getMass(PION), 2)) *
                                    (pow(deltaRestEnergy, 2) - pow(Particle::getMass(PROTON) - Particle::getMass(PION), 2))) /
                               (2 * deltaRestEnergy);

    // Setting the momentum in opposite directions
    proton.setRandomMomentumDirection(momentumMagnitude, *this);
    pion.setMomentum(-proton.px, -proton.py, -proton.pz);

    // 2. Lorentz boost to transform the momenta and energy from the rest frame to the lab frame
    auto [vx, vy, vz] = deltaParticle.velocity();
    proton.lorentzBoost(vx, vy, vz);
    pion.lorentzBoost(vx, vy, vz);

    // Set the position of the new particles to the position of the decaying DELTA
    proton.x = deltaParticle.x;
    proton.y = deltaParticle.y;
    proton.z = deltaParticle.z;

    pion.x = deltaParticle.x;
    pion.y = deltaParticle.y;
    pion.z = deltaParticle.z;

    proton.birthStep = currentStep;
    pion.birthStep = currentStep;

    return {proton, pion};
}

//-----------------------------------------------------------------------------

void Simulation::processInteractions(int currentStep)
{
    std::vector<Particle> newParticles;
    std::unordered_set<size_t> indicesToRemoveSet; // 使用unordered_set而不是vector
    auto &particles = box.getParticles();

    for (size_t i = 0; i < particles.size(); ++i)
    {
        Particle &currentParticle = particles[i];

        // 检查粒子是否已经被标记为删除
        if (indicesToRemoveSet.find(i) != indicesToRemoveSet.end())
        {
            continue;
        }

        if (currentParticle.type == DELTA && decayOccurs(currentParticle))
        {
            auto [proton, pion] = createProtonAndPionFromDecay(currentParticle, currentStep);

            indicesToRemoveSet.insert(i);

            newParticles.push_back(proton);
            newParticles.push_back(pion);
        }
        else if (currentParticle.type == PROTON || currentParticle.type == PION)
        {
            auto neighbors = box.getNeighboringParticles(currentParticle.x, currentParticle.y, currentParticle.z);
            for (Particle *neighbor : neighbors)
            {
                size_t neighborIndex = neighbor - &particles[0];

                if (indicesToRemoveSet.find(neighborIndex) != indicesToRemoveSet.end())
                {
                    continue;
                }

                if (&currentParticle != neighbor &&
                    ((currentParticle.type == PROTON && neighbor->type == PION) || (currentParticle.type == PION && neighbor->type == PROTON)))
                {
                    if (&currentParticle < neighbor && collisionOccurs(currentParticle, *neighbor, currentStep))
                    {
                        Particle delta = createDeltaFromCollision(currentParticle, *neighbor, currentStep);

                        indicesToRemoveSet.insert(i);

                        indicesToRemoveSet.insert(neighborIndex);

                        newParticles.push_back(delta);
                        break;
                    }
                }
            }
        }
    }

    // 删除粒子
    std::vector<size_t> indicesToRemove(indicesToRemoveSet.begin(), indicesToRemoveSet.end());
    std::sort(indicesToRemove.rbegin(), indicesToRemove.rend());
    for (size_t index : indicesToRemove)
    {
        box.removeParticle(index);
    }

    // 添加新的粒子
    for (const Particle &p : newParticles)
    {
        box.addParticle(p);
    }

    indicesToRemoveSet.clear();
    indicesToRemove.clear();
    newParticles.clear();
}

//
void Simulation::run()
{

    outFile.open("particles.csv", std::ios::app);

    if (!outFile.is_open())
    {
        std::cerr << "Unable to open file for writing." << std::endl;
        return;
    }

    int steps = static_cast<int>(totalTime / dt);
    for (int event = 0; event < numEvents; ++event)
    {
        std::cout << "events:==" << event << std::endl;
        //        outputEventHeader(event);

        initializeParticles();
        for (int i = 0; i < steps; ++i)
        {
            box.assignParticlesToGrid();
            processInteractions(i);

            updateParticlePositions();
            outputCurrentState(i);
        }
    }
    outFile.close();
}

void Simulation::updateParticlePositions()
{
    for (Particle &particle : box.getParticles())
    {
        particle.move(dt);
    }
    box.applyPeriodicBoundaryConditions();
}

void Simulation::outputEventHeader(int eventNumber)
{

    if (outFile.is_open())
    {
        outFile << "Event Number: " << eventNumber << std::endl;
        outFile << "Parameters: " << std::endl;
        outFile << "Time Step: " << dt << ", Total Time: " << totalTime;
        outFile << ", Box Size: " << boxSize << ", Desired Cell Size: " << desiredCellSize;
        outFile << ", Protons: " << numProtons << ", Pions: " << numPions;
        outFile << ", kT: " << kT << ", Seed: " << seed << std::endl;
        outFile << "---------------------------------------" << std::endl;

        // Optionally print to console
        std::cout << "Starting Event Number " << eventNumber << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}

void Simulation::outputCurrentState(int step)
{
    std::map<ParticleType, int> particleCounts = box.countParticles();
    double totalEnergy = 0.0;

    // Calculate the total energy of all particles
    for (const Particle &particle : box.getParticles())
    {
        totalEnergy += particle.energy;
    }
    bool check = false;
    if (outFile.is_open())
    {
        // If it's the first time step, write the header
        if (step == 0)
        {
            outFile << "Time step,PROTON,PION,DELTA,Total Energy" << std::endl;
        }

        // Write data
        outFile << step << ","
                << particleCounts[PROTON] << ","
                << particleCounts[PION] << ","
                << particleCounts[DELTA] << ","
                << totalEnergy
                << std::endl;

        if (particleCounts[PROTON] + particleCounts[PION] + 2 * particleCounts[DELTA] != numPions + numProtons)
        {
            check = true;
        }
        if (check)
        {
            std::cerr << "Error: Particle conservation violated at step " << step << std::endl;
            exit(-1);
        }
    }
    else
    {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}
