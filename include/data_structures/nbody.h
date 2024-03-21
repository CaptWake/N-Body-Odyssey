#ifndef NBODY_H
#define NBODY_H

#include <vector>

#include "particle.h"

class NBody {
 public:
  double gravitationalConstant;
  // Constructors//
  NBody();
  NBody(const std::vector<Particle>& particles, double gravitationalConstant);

  // Getters//
  const std::vector<Particle>& GetParticles() const { return particles; }

  // Setters//
  void SetParticles(const std::vector<Particle>& particles) {
    this->particles = particles;
  }

  // Update//
  void update(double timestep);

  friend std::ostream& operator<<(std::ostream& os, const NBody& nbody);

  void addParticle(Particle particle);

 protected:
  std::vector<Particle> particles;

 private:
  vec3 CalculateGravitationalForce(Particle& particle1,
                                   Particle& particle2) const;
};

#endif
