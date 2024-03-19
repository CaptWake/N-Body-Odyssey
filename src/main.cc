#include "../data_structures/particle.h"
#include <vector>

vec3 CalculateGravitationalForce(const Particle& particle1, const Particle& particle2, double gravitational_constant) {
  // Gravitational constant (G)
  const double G = gravitational_constant;

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


void AllPairsSequential(std::vector<Particle>& particles, double timestep, double gravitational_constant) {
  auto num_particles = particles.size();

  // Loop through all particle pairs (avoid double counting)
  for (int i = 0; i < num_particles; ++i) {
    for (int j = i + 1; j < num_particles; ++j) {
      // Calculate gravitational force on particle i due to particle j
      vec3 force_on_i = CalculateGravitationalForce(particles[i], particles[j], gravitational_constant);

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

int main() {

  std::vector<Particle> particles;

  auto mercury = Particle(1.6601367952719304e7,
                          vec3(3.067925455119497e1, -2.695851754307215e-01, -5.017408038016835e-02),
                          vec3(1.300791785131080e-02, 2.245507195129482e-02, 6.414809034847098e-04),
                          {0, 0, 0}, 20, 5.43);

  std::cout << mercury << std::endl;

  auto venus = Particle(2.4478383396645447e-06,
                        vec3(7.248985943229426e-01,  2.250842837066930e-02, -4.152252796285734e-02),
                        vec3(-7.115663662455383e-04,  2.012442278126239e-02,  3.170051245998433e-04), {0, 0, 0}, 20, 5.24);

  std::cout << venus << std::endl;

 // std::cout << CalculateGravitationalForce(mercury, venus, 6.6743e-11).length() << std::endl;

  std::cout << "================= TEST n.1 WITH ONE DIMENSION =================" << std::endl;

  // Taking example from the following resource https://www.sciencefacts.net/gravitational-force.html
  auto earth = Particle(6e24, vec3(1.5e11, 0,0), vec3(0,0,0));
  auto sun = Particle(2e30, vec3(0,0,0), vec3(0,0,0));
  auto f_hat = CalculateGravitationalForce(earth, sun, 6.6743e-11).length();

  std::cout << "Estimated Gravitational Force: " << f_hat << std::endl;

  auto f_true = 3.5e22;
  std::cout << "Expected Gravitational Force: " << f_true << std::endl;

  std::cout << "================= TEST n.2 WITH ONE DIMENSION =================" << std::endl;

  earth = Particle(5.972e24, vec3(0, 0,0), vec3(0,0,0));
  auto moon = Particle(7.348e22, vec3( 3.844e8,0,0), vec3(0,0,0));
  f_hat = CalculateGravitationalForce(earth, moon, 6.6743e-11).length();

  std::cout << "Estimated Gravitational Force: " << f_hat << std::endl;

  std::cout << "================= TEST ALL PAIRS ALGORITHM =================" << std::endl;

  particles.emplace_back(6.66666667, vec3(0.27626589, -1.85462808, 0.62390111), vec3(-0.43778315, 2.171257, 1.15231025));
  particles.emplace_back(6.66666667, vec3(1.14531129, 1.03719047, 1.88663893), vec3(-1.81881234, -0.13804934, 0.53983961));
  particles.emplace_back(6.66666667, vec3(-0.11169829, -0.36210134, 0.14867505), vec3(-1.77528229, 1.31487654, -0.47344805));

  AllPairsSequential(particles, 0.01, 1.0);
  for (auto& particle : particles) {
    std::cout << "Position: " << particle.GetPosition() << std::endl;
    std::cout << "Velocity: " << particle.GetVelocity() << std::endl;
  }
  return 0;
}