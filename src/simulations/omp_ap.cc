#include "omp_ap.h"
#include <omp.h>
#include <cmath>
#include "nbody_helpers.h"
#include "time_utils.h"
#include <iostream>

// copyright NVIDIA
void OMPAPUpdate(const uint64_t n, float *m, float *p , float *v, const float dt) {

  #pragma omp parallel for schedule(runtime)
  for (uint64_t i = 0; i < n * 3; i += 3) {
    float fx = 0.0f;
    float fy = 0.0f;
    float fz = 0.0f;
    for (uint64_t j = 0; j < n * 3; j += 3) {
      auto m2_id = j / 3;
      // compute distance pair
      auto dx = p[j] - p[i];

      auto dy = p[j + 1] - p[i + 1];
      auto dz = p[j + 2] - p[i + 2];

      auto d = dx * dx + dy * dy + dz * dz + _SOFTENING*_SOFTENING;
      auto d_inv = 1.0f / sqrtf(d);
      auto d_inv3 = d_inv * d_inv * d_inv;

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
void OMPAPSimulate(uint64_t n, float dt, float tEnd, uint64_t seed){

  float *m = new float[n];
  float *p = new float[3 * n];
  float *v = new float[3 * n];
  float *a = new float[3 * n];

  // Init Bodies
  InitAos(n, m, p, v, a);

  // Simulation Loop
  for (float t = 0.0f; t < tEnd; t += dt){
    // Update Bodies
    OMPAPUpdate(n, m, p, v, dt);
    float ek = Ek(n, m, v);
    float ep = Ep(n, m, p);
    std::cout << "Etot: " <<ek+ep <<std::endl;
  }
}

int main (int argc, char **argv) {
  if (argc < 5) {
    std::cerr << "Must specify the number of bodies, schedule type, chunk size and number of threads" << std::endl;
    exit(1);
  }
  srand(0);
  SetScheduleType(argv[2], atoi(argv[3]));
  SetNumThread(atoi(argv[4]));
  TIMERSTART(simulation)
  OMPAPSimulate(std::stoul(argv[1]), 0.01, 0.1, 0);
  TIMERSTOP(simulation)
}
