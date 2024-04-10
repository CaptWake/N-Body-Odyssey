#ifndef OMP_AP_H_
#define OMP_AP_H_

#include <cstdint>
#include <vector>
#include <string>
#include <random>

#include "nbody.h"
#include "fileIO.h"


class OmpAP : NBody {
 public:
  OmpAP() = default;

  OmpAP(const std::string& fname, const int num_threads, const std::string& schedule_type, const int chunk_size) {
    std::vector<float> m, v, p;
    this->num_threads = num_threads;
    this->schedule_type = schedule_type;
    this->chunk_size = chunk_size;

    this->n_bodies = ReadCSVConfiguration(fname, m, p, v, this->G);

    this->masses = std::move(m);
    this->positions = std::move(p);
    this->velocities = std::move(v);
  }

  // generate random samples
  OmpAP(uint64_t n_bodies, const float grav_const, const int num_threads, const std::string& schedule_type, const int chunk_size) {
    static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
    // initialize Mersennes' twister using rd to generate the seed
    static std::mt19937 engine{0};//rd()};
    std::uniform_real_distribution<float> density(-1, 1);
    const uint64_t n_coords = n_bodies*3;

    this->n_bodies = n_bodies;
    this->G = grav_const;
    this->masses.reserve(n_bodies);
    this->positions.reserve(n_coords);
    this->velocities.reserve(n_coords);

    this->num_threads = num_threads;
    this->schedule_type = schedule_type;
    this->chunk_size = chunk_size;

    for (uint64_t i = 0; i < n_coords; i+=3) {
      this->masses.push_back(density(engine));

      this->positions.push_back(density(engine));
      this->positions.push_back(density(engine));
      this->positions.push_back(density(engine));

      this->velocities.push_back(density(engine));
      this->velocities.push_back(density(engine));
      this->velocities.push_back(density(engine));
    }
  }

  //Move Constructor
  OmpAP& operator=(OmpAP&& old) noexcept {
    n_bodies=old.n_bodies;
    masses=std::move(old.masses);
    velocities=std::move(old.velocities);
    positions=std::move(old.positions);
    G=old.G;
    num_threads=old.num_threads;
    schedule_type=old.schedule_type;
    chunk_size=old.chunk_size;
    return *this;
  }

  // Function to export bodies information to CSV file
  void LogsToCSV(const std::string& filename) const;

  // Update//
  void Update(float dt) override;

  friend std::ostream& operator<<(std::ostream& os, const OmpAP& nbody);

  std::vector<float> masses{};
  std::vector<float> positions{};
  std::vector<float> velocities{};
  uint64_t n_bodies{};
  std::string schedule_type;
  int chunk_size{}, num_threads{};
  float G{};

};

#endif
