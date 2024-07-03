#include "simulations/omp_ap.h"

#include <omp.h>

#include <cmath>
#include <iostream>

#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

/**
 * @brief Updates the positions and velocities of bodies using OpenMP
 * parallelization.
 *
 * @tparam T Floating point type (float or double).
 * @param n Total number of bodies.
 * @param m Array of masses.
 * @param p Array of positions.
 * @param v Array of velocities.
 * @param dt Time step for the simulation.
 */
template <typename T>
void OMPAPUpdate(const int n, T *m, T *p, T *v, const T dt) {
#pragma omp parallel for schedule(runtime)
  for (int i = 0; i < n * 3; i += 3) {
    T fx = 0.0f;
    T fy = 0.0f;
    T fz = 0.0f;
    for (int j = 0; j < n * 3; j += 3) {
      int m2_id = j / 3;
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

  for (int i = 0; i < n * 3; i += 3) {
    p[i] += v[i] * dt;
    p[i + 1] += v[i + 1] * dt;
    p[i + 2] += v[i + 2] * dt;
  }
}

/**
 * @brief Simulates the n-body problem using OpenMP.
 *
 * @tparam T Floating point type (float or double).
 * @param n Total number of bodies.
 * @param dt Time step for the simulation.
 * @param tEnd End time for the simulation.
 * @param seed Random seed for initialization.
 */
template <typename T>
void OMPAPSimulate(int n, T dt, T tEnd, int seed) {
  T *m = new T[n];
  T *p = new T[3 * n];
  T *v = new T[3 * n];

  // Init Bodies
  InitAos<T>(n, m, p, v, seed);

#ifdef MONITOR_ENERGY
  T ek = Ek<T>(n, m, v);
  T ep = Ep<T>(n, m, p);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif

  TIMERSTART(simulation)
  // Simulation Loop
  for (T t = 0.0f; t < tEnd; t += dt) {
    // Update Bodies
    OMPAPUpdate<T>(n, m, p, v, dt);

#ifdef MONITOR_ENERGY
    ek = Ek<T>(n, m, v);
    ep = Ep<T>(n, m, p);
    std::cout << "Etot: " << ek + ep << std::endl;
#endif
  }

#ifdef MONITOR_ENERGY
  ek = Ek<T>(n, m, v);
  ep = Ep<T>(n, m, p);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif

  TIMERSTOP(simulation)
  delete[] m;
  delete[] p;
  delete[] v;
}

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cerr << "Must specify the number of bodies, schedule type, chunk "
                 "size, number of threads and seed (optional)"
              << std::endl;
    exit(1);
  }
  int nbody = atoi(argv[1]);
  char *scheduleType = argv[2];
  int blockSize = atoi(argv[3]);
  int numThread = atoi(argv[4]);
  int seed = 0;

  if (argc == 6) seed = atoi(argv[5]);

  SetScheduleType(scheduleType, blockSize);
  SetNumThread(numThread);

  OMPAPSimulate<MY_T>(nbody, 0.01, 1, seed);
}
