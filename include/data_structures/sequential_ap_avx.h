#ifndef SEQUENTIAL_AP_AVX_H_
#define SEQUENTIAL_AP_AVX_H_

#include <cstdint>
#include <iostream>
#include <random>

#include "fileIO.h"
#include "nbody.h"

class SequentialAPAVX : NBody {
 public:
  SequentialAPAVX() {
    G = 1;
    m = nullptr;
    px = nullptr;
    py = nullptr;
    pz = nullptr;
    vx = nullptr;
    vy = nullptr;
    vz = nullptr;
    n_bodies = 0;
  }

  explicit SequentialAPAVX(const std::string &fname) {
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
  }

  // generate random samples
  SequentialAPAVX(uint64_t n_bodies, float grav_const) {
    static std::random_device rd;  // random device engine, usually based on
                                   // /dev/random on UNIX-like systems
    // initialize Mersennes' twister using rd to generate the seed
    static std::mt19937 engine{0};  // rd()};
    std::uniform_real_distribution<float> density(-1, 1);

    this->n_bodies = n_bodies;

    // rounds the number of slots to multiple of 8
    uint64_t n_slots = 8 * floor(n_bodies / 8 + 1);

    this->m = static_cast<float *>(_mm_malloc(n_slots * sizeof(float), 32));
    this->px = static_cast<float *>(_mm_malloc(n_slots * sizeof(float), 32));
    this->py = static_cast<float *>(_mm_malloc(n_slots * sizeof(float), 32));
    this->pz = static_cast<float *>(_mm_malloc(n_slots * sizeof(float), 32));
    this->vx = static_cast<float *>(_mm_malloc(n_slots * sizeof(float), 32));
    this->vy = static_cast<float *>(_mm_malloc(n_slots * sizeof(float), 32));
    this->vz = static_cast<float *>(_mm_malloc(n_slots * sizeof(float), 32));
    this->G = grav_const;

    for (uint64_t i = 0; i < n_bodies; ++i) {
      this->m[i] = density(engine);
      this->px[i] = density(engine);
      this->py[i] = density(engine);
      this->pz[i] = density(engine);
      this->vx[i] = density(engine);
      this->vy[i] = density(engine);
      this->vz[i] = density(engine);
    }
  }

  // Move Constructor
  SequentialAPAVX &operator=(SequentialAPAVX &&old) noexcept {
    n_bodies = old.n_bodies;
    m = old.m;
    px = old.px;
    py = old.py;
    pz = old.pz;
    vx = old.vx;
    vy = old.vy;
    vz = old.vz;
    G = old.G;
    old.m = nullptr;
    old.px = nullptr;
    old.py = nullptr;
    old.pz = nullptr;
    old.vx = nullptr;
    old.vy = nullptr;
    old.vz = nullptr;
    return *this;
  }

  // Function to export bodies information to CSV file
  void LogsToCSV(const std::string &filename) const;

  // Update//
  void Update(float dt) override;

  friend std::ostream &operator<<(std::ostream &os,
                                  const SequentialAPAVX &nbody);

  ~SequentialAPAVX() {
    free(m);
    free(px);
    free(py);
    free(pz);
    free(vx);
    free(vy);
    free(vz);
  }

  float *m, *px, *py, *pz, *vx, *vy, *vz, G;
  uint64_t n_bodies;
};

#endif
