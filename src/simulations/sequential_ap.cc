#include "simulations/sequential_ap.h"

#include <cmath>

#include "utilities/integrators.h"
#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

template <typename T>
void computeForce(int n, int i, const T *m, const T *p, T *a) {
  //         a[k] = -p[k] / (r2*sqrtr2);
  T *ai = a + 3 * i;
  ai[0] = 0.0f;
  ai[1] = 0.0f;
  ai[2] = 0.0f;
  T px = p[3 * i + 0];
  T py = p[3 * i + 1];
  T pz = p[3 * i + 2];
  T dx, dy, dz, D;
  int j;
  for (j = 0; j < i; ++j) {
    // really dx is other way around, be this way we can avoid -1.0* later.
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
  for (j = i + 1; j < n; ++j) {
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

// copyright NVIDIA
template <typename T>
void SequentialAPUpdate(const int n, T *m, T *p, T *v, const T dt) {
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

// Euler step https://en.wikipedia.org/wiki/File:Euler_leapfrog_comparison.gif//
template <typename T>
void SequentialAPSimulateV1(int n, T dt, T tEnd) {
  T *m = new T[n];
  T *p = new T[3 * n];
  T *v = new T[3 * n];

  // Init Bodies
  InitAos<T>(n, m, p, v);

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
#ifdef MONITOR_MOMENTUM
  std::array<T, 3> L = AngularMomentum<T>(n, m, p, v);
  T mod = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
  std::cout << "L: " << "(" << L[0] << ", " << L[1] << ", " << L[2] << ") mod: "<< mod << std::endl;
#endif
  }
#ifdef MONITOR_ENERGY
  ek = Ek<T>(n, m, v);
  ep = Ep<T>(n, m, p);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif
#ifdef MONITOR_MOMENTUM
  std::array<T, 3> L = AngularMomentum<T>(n, m, p, v);
  T mod = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
  std::cout << "L: " << "(" << L[0] << ", " << L[1] << ", " << L[2] << ") mod: "<< mod << std::endl;
#endif
  TIMERSTOP(simulation)
}


// Kick-drift-kick //
template <typename T>
void SequentialAPSimulateV2(int n, T dt, T tEnd) {
  T *m = new T[n];
  T *p = new T[3 * n];
  T *v = new T[3 * n];
  T *a = new T[3 * n];

  // Init Bodies
  InitAos<T>(n, m, p, v, a);
#ifdef MONITOR_ENERGY
  T ek = Ek<T>(n, m, v);
  T ep = Ep<T>(n, m, p);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif
#ifdef MONITOR_MOMENTUM
  std::array<T, 3> L = AngularMomentum<T>(n, m, p, v);
  T mod = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
  std::cout << "L: " << "(" << L[0] << ", " << L[1] << ", " << L[2] << ") mod: "<< mod << std::endl;
#endif

  // Simulation Loop
  TIMERSTART(simulation)
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
#ifdef MONITOR_MOMENTUM
  L = AngularMomentum<T>(n, m, p, v);
  mod = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
  std::cout << "L: " << "(" << L[0] << ", " << L[1] << ", " << L[2] << ") mod: "<< mod << std::endl;
#endif
  }
  TIMERSTOP(simulation)
#ifdef MONITOR_ENERGY
  ek = Ek<T>(n, m, v);
  ep = Ep<T>(n, m, p);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif
#ifdef MONITOR_MOMENTUM
  L = AngularMomentum<T>(n, m, p, v);
  mod = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
  std::cout << "L: " << "(" << L[0] << ", " << L[1] << ", " << L[2] << ") mod: "<< mod << std::endl;
#endif
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  if (argc == 3)
    srand(atoi(argv[2]));
  else
    srand(0);
  SequentialAPSimulateV1<MY_T>(atoi(argv[1]), 0.01, 10);
}