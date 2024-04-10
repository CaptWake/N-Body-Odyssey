#ifndef OMP_BH_H_
#define OMP_BH_H_

#include "nbody.h"
#include "vec3.h"
#include <utility>
#include <vector>
#include <random>

class OmpBH : NBody {
 public:
  OmpBH() = default;

  // generate random samples
  OmpBH(const uint64_t n_bodies, const float grav_const, const float theta, const int num_threads, const std::string& schedule_type="dynamic", const int chunk_size=1) {
    static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
    // initialize Mersennes' twister using rd to generate the seed
    static std::mt19937 engine{0};//rd()};
    std::uniform_real_distribution<float> density(-1, 1);

    this->G = grav_const;
    this->theta = theta;
    this->num_threads = num_threads;
    this->schedule_type = schedule_type;
    this->chunk_size = chunk_size;

    this->m.reserve(n_bodies);
    this->p.reserve(n_bodies);
    this->v.reserve(n_bodies);

    for (uint64_t i = 0; i < n_bodies; i++) {
      this->m.push_back(density(engine));
      this->p.emplace_back(density(engine), density(engine), density(engine));
      this->v.emplace_back(density(engine), density(engine), density(engine));
    }
  }

  void Update(float dt) override;
  void LogsToCSV(const std::string &filename);

  std::vector<vec3> p, v;
  std::vector<float> m;
  float G{}, theta{};
  std::string schedule_type;
  int chunk_size{}, num_threads{};
};

#endif
