#include "sequential_ap.h"

#include "data_structures/nbody.h"

// copyright NVIDIA

void SequentialAP::Update(const float dt) {

  std::cout << std::endl << std::endl;
  for (uint64_t i=0; i < this->n_bodies * 3; i+=3) {
    float fx = 0.0f;
    float fy = 0.0f;
    float fz = 0.0f;
    for (uint64_t j=0; j < this->n_bodies * 3; j+=3) {
      auto m2_id= j/3;
      // compute distance pair
      auto dx = this->positions[j] - this->positions[i];
      auto dy = this->positions[j+1] - this->positions[i+1];
      auto dz = this->positions[j+2] - this->positions[i+2];

      auto d = dx * dx + dy * dy + dz * dz + 1e-9f;
      auto d_inv = 1.0f / sqrtf(d);
      auto d_inv3 = d_inv * d_inv * d_inv;

      fx += d_inv3 * this->masses[m2_id] * dx;
      fy += d_inv3 * this->masses[m2_id] * dy;
      fz += d_inv3 * this->masses[m2_id] * dz;
    }

    this->velocities[i]   += fx * dt;
    this->velocities[i+1] += fy * dt;
    this->velocities[i+2] += fz * dt;

  }

  for (uint64_t i=0; i < this->n_bodies * 3; i+=3) {

    this->positions[i]   += this->velocities[i] * dt;
    this->positions[i+1] += this->velocities[i+1] * dt;
    this->positions[i+2] += this->velocities[i+2] * dt;
  }
}

std::ostream &operator<<(std::ostream &os, const SequentialAP &nbody) {
  std::cout << "Gravitational constant: " << nbody.G << std::endl;

  for (uint64_t i = 0; i < nbody.n_bodies; ++i) {
    std::cout << "Body " << i + 1 << ":\n";
    std::cout << "  Mass: " << nbody.masses[i] << std::endl;
    std::cout << "  Position (x, y, z): " << nbody.positions[i * 3] << ", " << nbody.positions[i * 3 + 1] << ", " << nbody.positions[i * 3 + 2] << std::endl;
    std::cout << "  Velocity (vx, vy, vz): " << nbody.velocities[i * 3] << ", " << nbody.velocities[i * 3 + 1] << ", " << nbody.velocities[i * 3 + 2] << std::endl;
  }
  return os;
}
