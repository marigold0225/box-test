#ifndef BOX_H
#define BOX_H

#include <vector>
#include "Particle.h"
#include <map>
#include <cstddef>
#include <cmath>

class Box
{
private:
    double size = 20.0;
    std::vector<Particle> particles;
    std::vector<std::vector<Particle *>> grid; // 添加格子数据结构

public:
    int gridSize;
    std::vector<Particle *> getNeighboringParticles(double x, double y, double z);
    std::map<ParticleType, int> countParticles() const;
    std::vector<Particle> &getParticles() { return particles; }
    void setSize(double s) { size = s; }
    double getSize() const { return size; }
    void addParticle(const Particle &p);
    void clearParticles();
    void applyPeriodicBoundaryConditions();
    void removeParticle(size_t index);

    void initializeGrid(double desiredCellSize);
    void assignParticlesToGrid();
    int getGridIndex(double coordinate) const;
    int getFlatIndex(int ix, int iy, int iz) const;
    int getGridSize() const { return gridSize; }

    double getCellVolume() const;
    const std::vector<std::vector<Particle *>> &getGrid() const { return grid; }
};

#endif // BOX_H
