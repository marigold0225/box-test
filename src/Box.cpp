#include "Box.h"

int Box::countParticles(ParticleType type) const
{
    int count = 0;
    for (const Particle &particle : particles)
    {
        if (particle.type == type)
        {
            count++;
        }
    }
    return count;
}

void Box::addParticle(const Particle &p)
{
    particles.push_back(p);
}

void Box::removeParticle(size_t index)
{
    if (index < particles.size())
    {
        particles.erase(particles.begin() + index);
    }
}

void Box::clearParticles()
{
    particles.clear();
}

void Box::applyPeriodicBoundaryConditions()
{
    for (Particle &p : particles)
    {
        if (p.x < 0)
            p.x += size;
        else if (p.x > size)
            p.x -= size;

        if (p.y < 0)
            p.y += size;
        else if (p.y > size)
            p.y -= size;

        if (p.z < 0)
            p.z += size;
        else if (p.z > size)
            p.z -= size;
    }
}

void Box::initializeGrid(double desiredCellSize)
{
    gridSize = static_cast<int>(std::ceil(size / desiredCellSize));
    grid.resize(gridSize * gridSize * gridSize);
}

void Box::assignParticlesToGrid()
{
    for (auto &cell : grid)
    {
        cell.clear();
    }
    for (Particle &particle : particles)
    {
        int ix = getGridIndex(particle.x);
        int iy = getGridIndex(particle.y);
        int iz = getGridIndex(particle.z);
        grid[getFlatIndex(ix, iy, iz)].push_back(&particle);
    }
}

int Box::getGridIndex(double coordinate) const
{
    int index = static_cast<int>(coordinate / size * gridSize);
    return std::min(index, gridSize - 1); // ensure the index is never equal to gridSize
}

int Box::getFlatIndex(int ix, int iy, int iz) const
{
    return ix * gridSize * gridSize + iy * gridSize + iz;
}

double Box::getCellVolume() const
{
    double cellSize = size / gridSize;
    return cellSize * cellSize * cellSize;
}

std::vector<Particle *> Box::getNeighboringParticles(double x, double y, double z)
{
    int ix = getGridIndex(x);
    int iy = getGridIndex(y);
    int iz = getGridIndex(z);

    std::vector<Particle *> neighboringParticles;

    for (int dx = -1; dx <= 1; ++dx)
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {
                int nx = ix + dx;
                int ny = iy + dy;
                int nz = iz + dz;

                if (nx < 0)
                    nx = gridSize - 1;
                else if (nx >= gridSize)
                    nx = 0;

                if (ny < 0)
                    ny = gridSize - 1;
                else if (ny >= gridSize)
                    ny = 0;

                if (nz < 0)
                    nz = gridSize - 1;
                else if (nz >= gridSize)
                    nz = 0;

                auto &cell = grid[getFlatIndex(nx, ny, nz)];
                neighboringParticles.insert(neighboringParticles.end(), cell.begin(), cell.end());
            }
        }
    }

    return neighboringParticles;
}
