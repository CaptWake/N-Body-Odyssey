#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_
#include "utilities/nbody_helpers.h"

static inline void performNBodyHalfStepA(uint64_t n, float dt,
                                         float* p,
                                         float* v,
                                         const float* a,
                                         const float* m){
  for (uint64_t i = 0; i < n; ++i) {
    //kick, drift
    v[3*i + 0] += 0.5f * a[3*i + 0] * dt;
    p[3*i + 0] += v[3*i + 0] * dt;
    v[3*i + 1] += 0.5f * a[3*i + 1] * dt;
    p[3*i + 1] += v[3*i + 1] * dt;
    v[3*i + 2] += 0.5f * a[3*i + 2] * dt;
    p[3*i + 2] += v[3*i + 2] * dt;
  }
}

static inline void performNBodyHalfStepB(uint64_t n, float dt,
                                         const float* p,
                                         float* v,
                                         const float* a,
                                         const float* m)
{
  for (uint64_t i = 0; i < n; ++i) {
    //kick
    v[3*i + 0] += 0.5f * a[3*i + 0] * dt;
    v[3*i + 1] += 0.5f * a[3*i + 1] * dt;
    v[3*i + 2] += 0.5f * a[3*i + 2] * dt;
  }
}

static inline void performNBodyStep(const int localN, const int n, float *m, float *p, float *p_, float *v, float *v_, const float dt) {
  for (int i = 0; i < localN*3; i += 3) {
    float fx = 0.0f;
    float fy = 0.0f;
    float fz = 0.0f;
    for (int j = 0; j < n * 3; j += 3) {
      auto m2_id = j / 3;
      // compute distance pair
      auto dx = p[j] - p_[i];
      auto dy = p[j + 1] - p_[i + 1];
      auto dz = p[j + 2] - p_[i + 2];

      auto d = dx * dx + dy * dy + dz * dz + _SOFTENING*_SOFTENING;
      auto d_inv = 1.0f / sqrtf(d);
      auto d_inv3 = d_inv * d_inv * d_inv;

      fx += d_inv3 * m[m2_id] * dx;
      fy += d_inv3 * m[m2_id] * dy;
      fz += d_inv3 * m[m2_id] * dz;
    }

    v_[i] += fx * dt;
    v_[i + 1] += fy * dt;
    v_[i + 2] += fz * dt;
  }

  for (int i = 0; i < localN*3; i += 3) {
    p_[i] += v_[i] * dt;
    p_[i + 1] += v_[i + 1] * dt;
    p_[i + 2] += v_[i + 2] * dt;
  }
}


#endif
