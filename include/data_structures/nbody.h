#ifndef NBODY_H
#define NBODY_H

#include "particle.h"
#include <vector>

class NBody {
 public:
  NBody(const std::vector<Particle> &particles, int end_time) : particles(particles), endTime(end_time) {}
  const std::vector<Particle> &GetParticles() const {
    return particles;
  }
  void SetParticles(const std::vector<Particle> &particles) {
    this->particles = particles;
  }
  int GetEndTime() const {
    return endTime;
  }
  void SetEndTime(int end_time) {
    endTime = end_time;
  }
 private:
  std::vector<Particle> particles;
  int endTime;
};

#endif
