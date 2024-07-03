#include "simulations/sequential_ap.h"

#include <cmath>
#include <iostream>

#include "utilities/integrators.h"
#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

/**
 * @brief Computes the force on a body using all other bodies.
 *
 * @tparam T Floating-point type (float or double).
 * @param n Total number of bodies.
 * @param i Index of the body for which the force is computed.
 * @param m Array of masses.
 * @param p Array of positions (AoS format).
 * @param a Array of accelerations (output).
 */
template <typename T>
void computeForce(int n, int i, const T *m, const T *p, T *a) {
  T *ai = a + 3 * i;
  ai[0] = 0.0f;
  ai[1] = 0.0f;
  ai[2] = 0.0f;
  T px = p[3 * i + 0];
  T py = p[3 * i + 1];
  T pz = p[3 * i + 2];
  T dx, dy, dz, D;

  for (int j = 0; j < i; ++j) {
    dx = p[3 * j] - px;
    dy = p[3 * j + 1] - py;
    dz = p[3 * j + 2] - pz;
    D = dx * dx + dy * dy + dz * dz;
    D += _SOFTENING * _SOFTENING;
    D = 1.0f / (D * sqrt(D));
    ai[0] += m[j] * dx * D;
    ai[1] += m[j] * dy * D;
    ai[2] += m[j] * dz * D;
  }

  for (int j = i + 1; j < n; ++j) {
    dx = p[3 * j] - px;
    dy = p[3 * j + 1] - py;
    dz = p[3 * j + 2] - pz;
    D = dx * dx + dy * dy + dz * dz;
    D += _SOFTENING * _SOFTENING;
    D = 1.0f / (D * sqrt(D));
    ai[0] += m[j] * dx * D;
    ai[1] += m[j] * dy * D;
    ai[2] += m[j] * dz * D;
  }
}

/**
 * @brief Updates positions and velocities of bodies using the Sequential
 * All-Pairs method.
 *
 * @tparam T Floating-point type (float or double).
 * @param n Total number of bodies.
 * @param m Array of masses.
 * @param p Array of positions (AoS format).
 * @param v Array of velocities (AoS format).
 * @param dt Time step for the simulation.
 */
template <typename T>
void SequentialAPUpdate(const int n, T *m, T *p, T *v, const T dt) {
  for (int i = 0; i < n * 3; i += 3) {
    T fx = 0.0f;
    T fy = 0.0f;
    T fz = 0.0f;
    for (int j = 0; j < n * 3; j += 3) {
      int m2_id = j / 3;
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
 * @brief Simulates the n-body problem using a basic Euler integration method.
 *
 * @tparam T Floating-point type (float or double).
 * @param n Total number of bodies.
 * @param dt Time step for the simulation.
 * @param tEnd End time for the simulation.
 * @param seed Seed for random number generator to initialize body states.
 */
template <typename T>
void SequentialAPSimulateV1(int n, T dt, T tEnd, int seed) {
  T *m = new T[n];
  T *p = new T[3 * n];
  T *v = new T[3 * n];

  // Initialize Bodies
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
    SequentialAPUpdate<T>(n, m, p, v, dt);
#ifdef MONITOR_ENERGY
    ek = Ek<T>(n, m, v);
    ep = Ep<T>(n, m, p);
    std::cout << "Etot: " << ek + ep << std::endl;
#endif
  }
  TIMERSTOP(simulation)
#ifdef MONITOR_ENERGY
  ek = Ek<T>(n, m, v);
  ep = Ep<T>(n, m, p);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif

  delete[] m;
  delete[] p;
  delete[] v;
}

/**
 * @brief Simulates the n-body problem using a kick-drift-kick integration
 * method.
 *
 * @tparam T Floating-point type (float or double).
 * @param n Total number of bodies.
 * @param dt Time step for the simulation.
 * @param tEnd End time for the simulation.
 */
template <typename T>
void SequentialAPSimulateV2(int n, T dt, T tEnd) {
  T *m = new T[n];
  T *p = new T[3 * n];
  T *v = new T[3 * n];
  T *a = new T[3 * n];

  // Initialize Bodies
  InitAos<T>(n, m, p, v, a);

#ifdef MONITOR_ENERGY
  T ek = Ek<T>(n, m, v);
  T ep = Ep<T>(n, m, p);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif

  TIMERSTART(simulation)
  // Simulation Loop
  for (T t = 0.0f; t < tEnd; t += dt) {
    // Update Bodies
    performNBodyHalfStepA<T>(n, dt, p, v, a, m);
    for (int i = 0; i < n; ++i) {
      computeForce<T>(n, i, m, p, a);
    }
    performNBodyHalfStepB<T>(n, dt, p, v, a, m);
#ifdef MONITOR_ENERGY
    ek = Ek<T>(n, m, v);
    ep = Ep<T>(n, m, p);
    std::cout << "Etot: " << ek + ep << std::endl;
#endif
  }
  TIMERSTOP(simulation)
#ifdef MONITOR_ENERGY
  ek = Ek<T>(n, m, v);
  ep = Ep<T>(n, m, p);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif

  delete[] m;
  delete[] p;
  delete[] v;
  delete[] a;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  int nbody = atoi(argv[1]);
  int seed = 0;
  if (argc == 3) seed = atoi(argv[2]);
#ifndef OMP
  srand(seed);
#endif
  SequentialAPSimulateV1<MY_T>(nbody, 0.01, 1, seed);
}
