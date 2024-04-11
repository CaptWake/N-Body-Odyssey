#ifndef SEQUENTIALAP_H_
#define SEQUENTIALAP_H_

#include <cstdint>
#include <iostream>
#include <random>

#include "fileIO.h"
#include "nbody.h"

class SequentialAP : public NBody {
 public:
  SequentialAP() = default;

  explicit SequentialAP(const std::string& fname) {
    std::vector<float> m, v, p;
    this->n_bodies = ReadCSVConfigurationAOS(fname, m, p, v, this->G);
    this->masses = std::move(m);
    this->positions = std::move(p);
    this->velocities = std::move(v);
  }

  // generate random samples
  SequentialAP(uint64_t n_bodies, const float grav_const) {
    static std::random_device rd;  // random device engine, usually based on
                                   // /dev/random on UNIX-like systems
    // initialize Mersennes' twister using rd to generate the seed
    static std::mt19937 engine{0};  // rd()};
    std::uniform_real_distribution<float> density(-1, 1);
    const uint64_t n_coords = n_bodies * 3;

    this->n_bodies = n_bodies;
    this->G = grav_const;
    this->masses.reserve(n_bodies);
    this->positions.reserve(n_coords);
    this->velocities.reserve(n_coords);

    for (uint64_t i = 0; i < n_coords; i += 3) {
      this->masses.push_back(density(engine));

      this->positions.push_back(density(engine));
      this->positions.push_back(density(engine));
      this->positions.push_back(density(engine));

      this->velocities.push_back(density(engine));
      this->velocities.push_back(density(engine));
      this->velocities.push_back(density(engine));
    }
  }

  // Move Constructor
  SequentialAP& operator=(SequentialAP&& old) noexcept {
    n_bodies = old.n_bodies;
    masses = std::move(old.masses);
    velocities = std::move(old.velocities);
    positions = std::move(old.positions);
    G = old.G;
    return *this;
  }

  // Function to export bodies information to CSV file
  void LogsToCSV(const std::string& filename) const;

  // Update//
  void Update(float dt) override;

  friend std::ostream& operator<<(std::ostream& os, const SequentialAP& nbody);

  std::vector<float> masses{};
  std::vector<float> positions{};
  std::vector<float> velocities{};
  uint64_t n_bodies{};
  float G{};
};

#endif
