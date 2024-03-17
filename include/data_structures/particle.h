#ifndef PARTICLE_H
#define PARTICLE_H

#include <ostream>
#include <iomanip>
#include "../utilities/vec3.h"

class Particle {
 private:
  double mass;
  double radius;
  double density;
  vec3 position;
  vec3 velocity;
  vec3 acceleration;

 public:
  Particle(double mass,
           const vec3 &position,
           const vec3 &velocity,
           const vec3 &acceleration={0,0,0},
           double radius=1,
           double density=1)
      : mass(mass),
        radius(radius),
        density(density),
        position(position),
        velocity(velocity),
        acceleration(acceleration) {}

  double GetMass() const {
    return mass;
  }
  void SetMass(double mass) {
    this->mass = mass;
  }
  double GetRadius() const {
    return radius;
  }
  void SetRadius(double radius) {
    this->radius = radius;
  }
  double GetDensity() const {
    return density;
  }
  void SetDensity(double density) {
    this->density = density;
  }
  const vec3 &GetPosition() const {
    return position;
  }
  void SetPosition(const vec3 &position) {
    this->position = position;
  }
  const vec3 &GetVelocity() const {
    return velocity;
  }
  void SetVelocity(const vec3 &velocity) {
    this->velocity = velocity;
  }
  const vec3 &GetAcceleration() const {
    return acceleration;
  }
  void SetAcceleration(const vec3 &acceleration) {
    this->acceleration = acceleration;
  }
  friend std::ostream& operator<<(std::ostream& os, const Particle& particle) {
    os << "Particle Information:\n";
    os << "  Mass:          " << std::setw(15) << std::setprecision(5) << particle.GetMass() << " solar masses\n";
    os << "  Radius:        " << std::setw(15) << std::setprecision(5) << particle.GetRadius() << " Hill radius\n";
    os << "  Density:       " << std::setw(15) << std::setprecision(5) << particle.GetDensity() << " g/cm^3\n";
    os << "  Position (x,y,z):" << std::setw(15) << particle.GetPosition()[0] << ", "
       << std::setw(15) << particle.GetPosition()[1] << ", "
       << std::setw(15) << particle.GetPosition()[2] << " AU\n";
    os << "  Velocity (x,y,z):" << std::setw(15) << particle.GetVelocity()[0] << ", "
       << std::setw(15) << particle.GetVelocity()[1] << ", "
       << std::setw(15) << particle.GetVelocity()[2] << " AU/days\n";
    os << "  Acceleration (x,y,z):" << std::setw(15) << particle.GetAcceleration()[0] << ", "
       << std::setw(15) << particle.GetAcceleration()[1] << ", "
       << std::setw(15) << particle.GetAcceleration()[2] << " [unknown]\n";
    return os;
  }

};

#endif
