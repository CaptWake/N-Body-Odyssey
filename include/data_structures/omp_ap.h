//
// Created by ste on 09/04/24.
//

#ifndef OMP_AP_H_
#define OMP_AP_H_

#include <cstdint>
#include <vector>
#include <string>
#include <random>

#include "nbody.h"


class OmpAP : NBody {
 public:
  OmpAP() = default;

  // generate random samples
  OmpAP(uint64_t n_bodies, float grav_const, const int num_threads, const std::string& schedule_type="dynamic", const int chunk_size=1) {
    static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
    // initialize Mersennes' twister using rd to generate the seed
    static std::mt19937 engine{0};//rd()};
    std::uniform_real_distribution<float> density(-1, 1);

    this->n_bodies   = n_bodies;
    this->masses     = new float[n_bodies];
    this->positions  = new float[n_bodies*3];
    this->velocities = new float[n_bodies*3];
    this->G = grav_const;
    this->num_threads=num_threads;
    this->chunk_size=chunk_size;
    this->schedule_type=schedule_type;

    for (uint64_t i = 0; i < n_bodies*3; i+=3) {
      this->masses[i/3]     = density(engine);

      this->positions[i]    = density(engine);
      this->positions[i+1]  = density(engine);
      this->positions[i+2]  = density(engine);

      this->velocities[i]   = density(engine);
      this->velocities[i+1] = density(engine);
      this->velocities[i+2] = density(engine);
    }
  }

  //Move Constructor
  OmpAP& operator=(OmpAP&& old) noexcept {
    n_bodies=old.n_bodies;
    masses=old.masses;
    velocities=old.velocities;
    positions=old.positions;
    G=old.G;
    num_threads=old.num_threads;
    schedule_type=old.schedule_type;
    chunk_size=old.chunk_size;
    old.masses = nullptr;
    old.positions = nullptr;
    old.velocities = nullptr;
    return *this;
  }

  // Function to export bodies information to CSV file
  void LogsToCSV(const std::string& filename) const;

  // Update//
  void Update(float dt) override;

  friend std::ostream& operator<<(std::ostream& os, const OmpAP& nbody);

  ~OmpAP() {
    free(masses);
    free(positions);
    free(velocities);
  }

  float *masses{};
  float *positions{};
  float *velocities{};
  uint64_t n_bodies{};
  std::string schedule_type;
  int chunk_size{}, num_threads{};
  float G{};

};

#endif
