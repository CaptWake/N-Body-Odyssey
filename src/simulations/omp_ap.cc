#include "simulations/omp_ap.h"

#include <omp.h>

#include <cmath>
#include <iostream>

#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

// copyright NVIDIA
template<typename T>
void OMPAPUpdate(const uint64_t n, T *m, T *p, T *v,
                 const T dt) {
#pragma omp parallel for schedule(runtime)
  for (uint64_t i = 0; i < n * 3; i += 3) {
    T fx = 0.0f;
    T fy = 0.0f;
    T fz = 0.0f;
    for (uint64_t j = 0; j < n * 3; j += 3) {
      uint64_t m2_id = j / 3;
      // compute distance pair
      T dx = p[j] - p[i];

      T dy = p[j + 1] - p[i + 1];
      T dz = p[j + 2] - p[i + 2];

      T d = dx * dx + dy * dy + dz * dz + _SOFTENING * _SOFTENING;
      T d_inv = 1.0f / sqrt(d);
      T d_inv3 = d_inv * d_inv * d_inv;

      fx += d_inv3 * m[m2_id] * dx;
      fy += d_inv3 * m[m2_id] * dy;
      fz += d_inv3 * m[m2_id] * dz;
    }

    v[i] += fx * dt;
    v[i + 1] += fy * dt;
    v[i + 2] += fz * dt;
  }

  for (uint64_t i = 0; i < n * 3; i += 3) {
    p[i] += v[i] * dt;
    p[i + 1] += v[i + 1] * dt;
    p[i + 2] += v[i + 2] * dt;
  }
}

// Euler step https://en.wikipedia.org/wiki/File:Euler_leapfrog_comparison.gif//
template<typename T>
void OMPAPSimulate(uint64_t n, T dt, T tEnd, uint64_t seed) {
  T *m = new T[n];
  T *p = new T[3 * n];
  T *v = new T[3 * n];
  T *a = new T[3 * n];

  // Init Bodies
  InitAos<T>(n, m, p, v, a);

  // Simulation Loop
  for (T t = 0.0f; t < tEnd; t += dt) {
    // Update Bodies
    OMPAPUpdate<T>(n, m, p, v, dt);
    T ek = Ek<T>(n, m, v);
    T ep = Ep<T>(n, m, p);
    std::cout << "Etot: " << ek + ep << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cerr << "Must specify the number of bodies, schedule type, chunk size "
                 "and number of threads"
              << std::endl;
    exit(1);
  }
  srand(0);
  SetScheduleType(argv[2], atoi(argv[3]));
  SetNumThread(atoi(argv[4]));
  TIMERSTART(simulation)
  OMPAPSimulate<MY_T>(std::stoul(argv[1]), 0.01, 0.1, 0);
  TIMERSTOP(simulation)
}
