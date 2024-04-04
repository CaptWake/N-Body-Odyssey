#include "sequential_ap.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>

// copyright NVIDIA
void SequentialAP::Update(const float dt) {
  for (uint64_t i=0; i < this->n_bodies * 3; i+=3) {
    float fx = 0.0f;
    float fy = 0.0f;
    float fz = 0.0f;
    for (uint64_t j=0; j < this->n_bodies * 3; j+=3) {
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
void SequentialAP::LogsToCSV(const std::string &filename) const {
  std::ofstream file(filename, std::ios_base::app);
  if (file.is_open()) {
    uint64_t i;
    for (i = 0; i < n_bodies*3; i+=3) {
      file << this->positions[i] << "," << this->positions[i+1] << "," << this->positions[i+2];
      if (i < this->n_bodies * 3 - 3)
        file << ",";
    }
    file << std::endl;
    file.close();
  } else {
    std::cerr << "Unable to open file: " << filename << std::endl;
  }
}

void SequentialAP::LoadFromCSVConfiguration(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Cannot find configuration file at: " << filename << std::endl;
    exit(1);
  }

  std::string line;
  uint64_t i = 0;

  // Assuming the header is not present

  std::getline(file, line);
  std::istringstream iss(line);
  std::string token;
  std::getline(iss, token);
  // Get the gravitational constant (first line)
  this->G = std::stof(token);

  std::vector<float> m, p, v;
  while (std::getline(file, line)) {
    iss = std::istringstream(line);
    // Read particle data from CSV fields
    float mass, x, y, z, vx, vy, vz;
    // Read and parse comma-separated values
    std::getline(iss, token, ',');
    mass = std::stof(token);

    std::getline(iss, token, ',');
    x = std::stof(token);

    std::getline(iss, token, ',');
    y = std::stof(token);

    std::getline(iss, token, ',');
    z = std::stof(token);

    std::getline(iss, token, ',');
    vx = std::stof(token);

    std::getline(iss, token, ',');
    vy = std::stof(token);

    std::getline(iss, token, ',');
    vz = std::stof(token);

    m.push_back(mass);
    p.push_back(x);
    p.push_back(y);
    p.push_back(z);
    v.push_back(vx);
    v.push_back(vy);
    v.push_back(vz);
    i++;
  }
  file.close();

  // Update dynamic arrays
  this->masses = (float *)malloc(i * sizeof(float));
  this->positions = (float *)malloc(3*i * sizeof(float));
  this->velocities = (float *) malloc(3*i * sizeof(float));

  for (uint64_t k = 0; k < 3*i; ++k) {
    this->masses[k/3] = m[k/3];
    this->positions[k] = p[k];
    this->velocities[k] = v[k];
  }
  this->n_bodies = i;
}
