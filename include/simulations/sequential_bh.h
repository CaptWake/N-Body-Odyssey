#ifndef SEQUENTIAL_BH_H_
#define SEQUENTIAL_BH_H_

#include <vector>
#include "utilities/fileIO.h"
#include "utilities/vec3.h"
#include "utilities/nbody_helpers.h"

class SequentialBH {
 public:
  SequentialBH() = default;

  // generate random samples
  SequentialBH(const uint64_t n_bodies, const float theta) {
    this->theta = theta;

    float *m = new float[n_bodies];
    float *p = new float[3*n_bodies];
    float *v = new float[3*n_bodies];

    InitAos<float>(n_bodies, m, p, v);

    for (uint64_t i = 0; i < 3*n_bodies; i+=3) {
      this->m.push_back(m[i/3]);
      this->p.emplace_back(p[i], p[i+1], p[i+2]);
      this->v.emplace_back(v[i], v[i+1], v[i+2]);
    }

    delete[] m;
    delete[] p;
    delete[] v;
  }

  void Update(float dt);

  std::vector<vec3> p, v;
  std::vector<float> m;
  float theta{};
};

#endif
