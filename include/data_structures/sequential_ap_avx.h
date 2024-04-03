//
// Created by ste on 31/03/24.
//

#ifndef SEQUENTIAL_AP_AVX_H_
#define SEQUENTIAL_AP_AVX_H_


#include <cstdint>
#include <iostream>

#include "nbody.h"

class SequentialAPAVX: NBody {
 public:
  SequentialAPAVX() {
    G=1;
    m=nullptr;
    px=nullptr;
    py=nullptr;
    pz=nullptr;
    vx=nullptr;
    vy=nullptr;
    vz=nullptr;
    n_bodies=0;
  }

  explicit SequentialAPAVX(const std::string& fname) {
    LoadFromCSVConfiguration(fname);
  }

  SequentialAPAVX(uint64_t n_bodies, float grav_const) {
    this->n_bodies = n_bodies;
    this->m = new float[n_bodies];
    this->px = new float[n_bodies];
    this->py = new float[n_bodies];
    this->pz = new float[n_bodies];
    this->vx = new float[n_bodies];
    this->vy = new float[n_bodies];
    this->vz = new float[n_bodies];
    this->G = grav_const;
  }

  SequentialAPAVX(uint64_t n_bodies, float* m, float* px, float* py, float* pz, float* vx, float* vy, float* vz, float grav_const) {
    this->n_bodies = n_bodies;
    this->m = m;
    this->px = px;
    this->py = py;
    this->pz = pz;
    this->vx = vx;
    this->vy = vy;
    this->vz = vz;
    this->G = grav_const;
  }

  //Move Constructor
  SequentialAPAVX& operator=(SequentialAPAVX&& old) noexcept {
    n_bodies=old.n_bodies;
    m=old.m;
    px=old.px;
    py=old.py;
    pz=old.pz;
    vx=old.vx;
    vy=old.vy;
    vz=old.vz;
    G=old.G;
    old.m = nullptr;
    old.px = nullptr;
    old.py = nullptr;
    old.pz = nullptr;
    old.vx = nullptr;
    old.vy = nullptr;
    old.vz = nullptr;
    return *this;
  }

  void LoadFromCSVConfiguration(const std::string& filename);

  // Function to export bodies information to CSV file
  void LogsToCSV(const std::string& filename);

  // Update//
  void Update(float dt) override;

  friend std::ostream& operator<<(std::ostream& os, const SequentialAPAVX& nbody);

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
