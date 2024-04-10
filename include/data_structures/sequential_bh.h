#ifndef SEQUENTIAL_BH_H_
#define SEQUENTIAL_BH_H_

#include "nbody.h"
#include "vec3.h"
#include <vector>
#include <random>

class SequentialBH : NBody {
 public:
  SequentialBH() { this->G = 1; this->theta = 0.9f; }
  SequentialBH(const std::string &filename) {
    LoadFromCSVConfiguration(filename);
  }

  // generate random samples
  SequentialBH(const uint64_t n_bodies, const float grav_const, float theta) {
    static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
    // initialize Mersennes' twister using rd to generate the seed
    static std::mt19937 engine{0};//rd()};
    std::uniform_real_distribution<float> density(-1, 1);

    this->G = grav_const;
    this->theta = theta;

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
  void LoadFromCSVConfiguration(const std::string &filename);
  void LogsToCSV(const std::string &filename);

  std::vector<vec3> p, v;
  std::vector<float> m;
  float G, theta;
};

#endif
