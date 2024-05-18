#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_
#include "utilities/nbody_helpers.h"

template<typename T>
static inline void performNBodyHalfStepA(uint64_t n, T dt, T* p,
                                         T* v, const T* a,
                                         const T* m) {
  for (uint64_t i = 0; i < n; ++i) {
    // kick, drift
    v[3 * i + 0] += 0.5 * a[3 * i + 0] * dt;
    p[3 * i + 0] += v[3 * i + 0] * dt;
    v[3 * i + 1] += 0.5 * a[3 * i + 1] * dt;
    p[3 * i + 1] += v[3 * i + 1] * dt;
    v[3 * i + 2] += 0.5 * a[3 * i + 2] * dt;
    p[3 * i + 2] += v[3 * i + 2] * dt;
  }
}

template<typename T>
static inline void performNBodyHalfStepB(uint64_t n, T dt, const T* p,
                                         T* v, const T* a,
                                         const T* m) {
  for (uint64_t i = 0; i < n; ++i) {
    // kick
    v[3 * i + 0] += 0.5 * a[3 * i + 0] * dt;
    v[3 * i + 1] += 0.5 * a[3 * i + 1] * dt;
    v[3 * i + 2] += 0.5 * a[3 * i + 2] * dt;
  }
}

template<typename T>
static inline void performNBodyStep(const int localN, const int n, T* m,
                                    T* p, T* p_, T* v, T* v_,
                                    const T dt) {
  for (int i = 0; i < localN * 3; i += 3) {
    T fx = 0.0;
    T fy = 0.0;
    T fz = 0.0;
    for (int j = 0; j < n * 3; j += 3) {
      int m2_id = j / 3;
      // compute distance pair
      T dx = p[j] - p_[i];
      T dy = p[j + 1] - p_[i + 1];
      T dz = p[j + 2] - p_[i + 2];

      T d = dx * dx + dy * dy + dz * dz + _SOFTENING * _SOFTENING;
      T d_inv = 1.0f / sqrt(d);
      T d_inv3 = d_inv * d_inv * d_inv;

      fx += d_inv3 * m[m2_id] * dx;
      fy += d_inv3 * m[m2_id] * dy;
      fz += d_inv3 * m[m2_id] * dz;
    }

    v_[i] += fx * dt;
    v_[i + 1] += fy * dt;
    v_[i + 2] += fz * dt;
  }

  for (int i = 0; i < localN * 3; i += 3) {
    p_[i] += v_[i] * dt;
    p_[i + 1] += v_[i + 1] * dt;
    p_[i + 2] += v_[i + 2] * dt;
  }
}

#endif
