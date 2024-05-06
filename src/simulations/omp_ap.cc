#include "omp_ap.h"

#include <omp.h>

#include <cmath>
#include <fstream>
#include <sstream>

#include "sequential_ap.h"

void OmpAP::Update(const float dt) {
  if (this->schedule_type == "static") {
    omp_set_schedule(omp_sched_static, this->chunk_size);
  } else {
    omp_set_schedule(omp_sched_dynamic, this->chunk_size);
  }
  omp_set_num_threads(this->num_threads);

#pragma omp parallel for schedule(runtime)
  for (uint64_t i = 0; i < this->n_bodies * 3; i += 3) {
    float fx = 0.0f;
    float fy = 0.0f;
    float fz = 0.0f;

    // #pragma omp parallel for reduction (+: fx, fy, fz)
    for (uint64_t j = 0; j < this->n_bodies * 3; j += 3) {
      auto m2_id = j / 3;
      // compute distance pair

      auto dx = this->positions[j] - this->positions[i];

      auto dy = this->positions[j + 1] - this->positions[i + 1];
      auto dz = this->positions[j + 2] - this->positions[i + 2];

      auto d = dx * dx + dy * dy + dz * dz + 1e-9f;
      auto d_inv = 1.0f / sqrtf(d);
      auto d_inv3 = d_inv * d_inv * d_inv;

      fx += d_inv3 * this->masses[m2_id] * dx;
      fy += d_inv3 * this->masses[m2_id] * dy;
      fz += d_inv3 * this->masses[m2_id] * dz;
    }

    this->velocities[i] += fx * dt;
    this->velocities[i + 1] += fy * dt;
    this->velocities[i + 2] += fz * dt;
  }

  for (uint64_t i = 0; i < this->n_bodies * 3; i += 3) {
    this->positions[i] += this->velocities[i] * dt;
    this->positions[i + 1] += this->velocities[i + 1] * dt;
    this->positions[i + 2] += this->velocities[i + 2] * dt;
  }
}

std::ostream &operator<<(std::ostream &os, const OmpAP &nbody) {
  std::cout << "Gravitational constant: " << nbody.G << std::endl;

  for (uint64_t i = 0; i < nbody.n_bodies; ++i) {
    std::cout << "Body " << i + 1 << ":\n";
    std::cout << "  Mass: " << nbody.masses[i] << std::endl;
    std::cout << "  Position (x, y, z): " << nbody.positions[i * 3] << ", "
              << nbody.positions[i * 3 + 1] << ", "
              << nbody.positions[i * 3 + 2] << std::endl;
    std::cout << "  Velocity (vx, vy, vz): " << nbody.velocities[i * 3] << ", "
              << nbody.velocities[i * 3 + 1] << ", "
              << nbody.velocities[i * 3 + 2] << std::endl;
  }
  return os;
}

void OmpAP::LogsToCSV(const std::string &filename) const {
  std::ofstream file(filename, std::ios_base::app);
  if (file.is_open()) {
    uint64_t i;
    for (i = 0; i < n_bodies * 3; i += 3) {
      file << this->positions[i] << "," << this->positions[i + 1] << ","
           << this->positions[i + 2];
      if (i < this->n_bodies * 3 - 3) file << ",";
    }
    file << std::endl;
    file.close();
  } else {
    std::cerr << "Unable to open file: " << filename << std::endl;
  }
}
