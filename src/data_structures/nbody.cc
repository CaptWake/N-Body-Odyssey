#include "data_structures/nbody.h"

// PROBABLY WE CAN SET IT TO PRIVATE
// AND EXPOSE ONLY THE START METHOD ??
// INSTEAD OF USING THE SIMULATION CLASS?
void SequentialAPNBody::update(const double dt) {
  auto num_particles = bodies.size();

  // Loop through all particle pairs (avoid double counting)
  for (int i = 0; i < num_particles; ++i) {
    for (int j = 0; j < num_particles; ++j) {
      if (i != j) {
        // Calculate gravitational force on particle i due to particle j
        // vec3 force_on_i = CalculateGravitationalForce(bodies[i], bodies[j]);

        vec3 distance = bodies[j].GetPosition() - bodies[i].GetPosition();
        double distance_magnitude = distance.length();

        std::cout << distance << " " << sqrt(distance.x() * distance.x() + distance.y() * distance.y() + distance.z() * distance.z()) << std::endl;
        // Update acceleration of particle i
        auto acc = gravitationalConstant * bodies[j].GetMass() / (distance.squared_length()) * unit_vector(distance);
        bodies[i].SetAcceleration(acc);
      }
    }
  }

  // Update velocities and positions using timestep
  for (int i = 0; i < num_particles; ++i) {
    bodies[i].SetVelocity(bodies[i].GetVelocity() +
        dt * bodies[i].GetAcceleration() * 0.5);
    bodies[i].SetPosition(bodies[i].GetPosition() +
        dt * bodies[i].GetVelocity());

    //std::cout << "PARTICLE "<< i << std::endl;
    //std::cout << bodies[i] << std::endl;
  }

}

vec3 SequentialAPNBody::CalculateGravitationalForce(Body& particle1,
                                        Body& particle2) const {
  // Gravitational constant (G)
  const double G = gravitationalConstant;

  // Distance vector between bodies
  vec3 distance = particle2.GetPosition() - particle1.GetPosition();

  // Check for close encounters (avoid division by zero)
  double distance_magnitude = distance.length();

  // if (distance_magnitude < some_small_value) {
  //   distance = distance.normalize() * some_small_value; // Set minimum
  //   distance
  // }
  //  Gravitational force calculation
  return G * particle1.GetMass() * particle2.GetMass() /
         (distance_magnitude * distance_magnitude * distance_magnitude) *
         distance;
}
std::ostream& operator<<(std::ostream& os, const SequentialAPNBody& nbody) {
  os << "NBody bodies:\n";
  for (const Body& particle : nbody.GetBodies()) {
    os << particle << "\n";
  }
  return os;
}

void SequentialAPNBody::addBody(Body body) {
  this->bodies.push_back(body);
}
