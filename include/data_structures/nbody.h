#ifndef NBODY_H
#define NBODY_H

#include "particle.h"
#include <vector>

class NBody {
 public:
  //Constructors//
  NBody(const std::vector<Particle> &particles, double gravitationalConstant);

  //Getters//
  const std::vector<Particle> &GetParticles() const {return particles;}

  //Setters//
  void SetParticles(const std::vector<Particle> &particles) {this->particles = particles;}

  //Update//
  void update(double timestep);

  friend std::ostream& operator<<(std::ostream& os, const NBody& nbody);

  void addParticle(Particle particle);

 private:
  std::vector<Particle> particles;
  double gravitationalConstant;

  vec3 CalculateGravitationalForce(Particle& particle1,
                                     Particle& particle2) const;
};

#endif
