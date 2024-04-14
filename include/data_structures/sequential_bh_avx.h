#ifndef SEQUENTIAL_BH_AVX_H_
#define SEQUENTIAL_BH_AVX_H_

#include <random>
#include <vector>

#include "fileIO.h"
#include "nbody.h"
#include "vec3.h"

class SequentialBHAVX : NBody {
 public:
  SequentialBHAVX() = default;

  SequentialBHAVX(const std::string &fname, const float theta) {
    std::vector<float> _m, _px, _py, _pz, _vx, _vy, _vz;
    this->n_bodies = ReadCSVConfigurationSOA(
        fname, _m, _px, _py, _pz, _vx, _vy, _vz, this->G);
    vector_to_arr(_m, &this->m, true);
    vector_to_arr(_px, &this->px, true);
    vector_to_arr(_py, &this->py, true);
    vector_to_arr(_pz, &this->pz, true);
    vector_to_arr(_vx, &this->vx, true);
    vector_to_arr(_vy, &this->vy, true);
    vector_to_arr(_vz, &this->vz, true);
    this->theta = theta;
  }

  // generate random samples
  SequentialBHAVX(const uint64_t n_bodies, const float grav_const, float theta) {
    static std::random_device rd;  // random device engine, usually based on
                                   // /dev/random on UNIX-like systems
    // initialize Mersennes' twister using rd to generate the seed
    static std::mt19937 engine{0};  // rd()};
    std::uniform_real_distribution<float> density(-1, 1);

    this->G = grav_const;
    this->theta = theta;

    /*
    this->m.reserve(n_bodies);
    this->p.reserve(n_bodies);
    this->v.reserve(n_bodies);

    for (uint64_t i = 0; i < n_bodies; i++) {
      this->m.push_back(density(engine));
      this->p.emplace_back(density(engine), density(engine), density(engine));
      this->v.emplace_back(density(engine), density(engine), density(engine));
    }
    */

  }

  void Update(float dt) override;
  void LogsToCSV(const std::string &filename);

  float *px, *py, *pz, *vx, *vy, *vz, *m;
  float G{}, theta{};
  uint64_t n_bodies{};
};

#endif
