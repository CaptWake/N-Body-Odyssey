#include "data_structures/nbody.h"

NBody::NBody(const std::vector<Particle>& particles,
             double gravitationalConstant)
    : particles{ particles }
    , gravitationalConstant{ gravitationalConstant }{}

NBody::NBody()
    : particles{}
    , gravitationalConstant{1}{}

void NBody::update(const double timestep){
  auto num_particles = particles.size();

  // Loop through all particle pairs (avoid double counting)
  for (int i = 0; i < num_particles; ++i) {
    for (int j = i + 1; j < num_particles; ++j) {
      // Calculate gravitational force on particle i due to particle j
      vec3 force_on_i = CalculateGravitationalForce(particles[i], particles[j]);

      // Update acceleration of particle i
      particles[i].SetAcceleration(particles[i].GetAcceleration() + force_on_i / particles[i].GetMass());
    }
  }

  // Update velocities and positions using timestep
  for (int i = 0; i < num_particles; ++i) {
    particles[i].SetVelocity(particles[i].GetVelocity() + timestep * particles[i].GetAcceleration());
    particles[i].SetPosition(particles[i].GetPosition() + timestep * particles[i].GetVelocity());
  }
}

vec3 NBody::CalculateGravitationalForce(Particle& particle1,
                                          Particle& particle2) const{
  // Gravitational constant (G)
  const double G = gravitationalConstant;

  // Distance vector between particles
  vec3 distance = particle2.GetPosition() - particle1.GetPosition();

  // Check for close encounters (avoid division by zero)
  double distance_magnitude = distance.length();

  //if (distance_magnitude < some_small_value) {
  //  distance = distance.normalize() * some_small_value; // Set minimum distance
  //}
  // Gravitational force calculation
  return G * particle1.GetMass() * particle2.GetMass() / (distance_magnitude * distance_magnitude * distance_magnitude) * distance;
}
std::ostream& operator<<(std::ostream& os, const NBody& nbody) {
  os << "NBody particles:\n";
  for (const Particle& particle : nbody.particles) {
    os << particle << "\n";
  }
  return os;
}

void NBody::addParticle(Particle particle) {this->particles.push_back(particle);}
