#include "simulations/sequential_ap_avx.h"

#include <iostream>
#include "utilities/avx.h"
#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

// copyright NVIDIA

#ifdef FLOAT
/**
 * @brief Update the positions and velocities of bodies using AVX for float precision.
 * 
 * @param n Number of bodies.
 * @param m Array of masses.
 * @param px Array of x positions.
 * @param py Array of y positions.
 * @param pz Array of z positions.
 * @param vx Array of x velocities.
 * @param vy Array of y velocities.
 * @param vz Array of z velocities.
 * @param dt Time step for the simulation.
 */
void SequentialAPAVXUpdate(const int n, float *m, float *px, float *py,
                           float *pz, float *vx, float *vy, float *vz,
                           const float dt) {
  static __m256 DT = _mm256_set1_ps(dt);
  static __m256 s = _mm256_set1_ps(_SOFTENING);

  for (int i = 0; i < n; ++i) {
    __m256 Fx = _mm256_set1_ps(0.0f);
    __m256 Fy = _mm256_set1_ps(0.0f);
    __m256 Fz = _mm256_set1_ps(0.0f);

    for (int j = 0; j < n; j += 8) {  // exploit SIMD computing blocks of 8 pairs each time
      const __m256 Xi = _mm256_broadcast_ss(px + i);
      const __m256 Yi = _mm256_broadcast_ss(py + i);
      const __m256 Zi = _mm256_broadcast_ss(pz + i);

      const __m256 Xj = _mm256_loadu_ps(px + j);
      const __m256 Yj = _mm256_loadu_ps(py + j);
      const __m256 Zj = _mm256_loadu_ps(pz + j);

      const __m256 X = _mm256_sub_ps(Xj, Xi);
      const __m256 Y = _mm256_sub_ps(Yj, Yi);
      const __m256 Z = _mm256_sub_ps(Zj, Zi);

      __m256 M = _mm256_loadu_ps(m + j);

      const __m256 D = _mm256_fmadd_ps(X, X, _mm256_fmadd_ps(Y, Y, _mm256_fmadd_ps(Z, Z, _mm256_mul_ps(s, s))));
      const __m256 D_inv = _mm256_rsqrt_ps(D);
      const __m256 D_inv3 = _mm256_mul_ps(D_inv, _mm256_mul_ps(D_inv, D_inv));

      Fx = _mm256_fmadd_ps(D_inv3, _mm256_mul_ps(M, X), Fx);
      Fy = _mm256_fmadd_ps(D_inv3, _mm256_mul_ps(M, Y), Fy);
      Fz = _mm256_fmadd_ps(D_inv3, _mm256_mul_ps(M, Z), Fz);
    }
    const __m256 Vx = _mm256_mul_ps(Fx, DT);
    const __m256 Vy = _mm256_mul_ps(Fy, DT);
    const __m256 Vz = _mm256_mul_ps(Fz, DT);
    vx[i] += hsum_avx(Vx);
    vy[i] += hsum_avx(Vy);
    vz[i] += hsum_avx(Vz);
  }

  for (int i = 0; i < n; i += 8) {
    const __m256 X = _mm256_loadu_ps(px + i);
    const __m256 Y = _mm256_loadu_ps(py + i);
    const __m256 Z = _mm256_loadu_ps(pz + i);

    const __m256 Vx = _mm256_loadu_ps(vx + i);
    const __m256 Vy = _mm256_loadu_ps(vy + i);
    const __m256 Vz = _mm256_loadu_ps(vz + i);

    const __m256 Xn = _mm256_fmadd_ps(Vx, DT, X);
    const __m256 Yn = _mm256_fmadd_ps(Vy, DT, Y);
    const __m256 Zn = _mm256_fmadd_ps(Vz, DT, Z);

    _mm256_store_ps(px + i, Xn);
    _mm256_store_ps(py + i, Yn);
    _mm256_store_ps(pz + i, Zn);
  }
}
#else
/**
 * @brief Update the positions and velocities of bodies using AVX for double precision.
 * 
 * @param n Number of bodies.
 * @param m Array of masses.
 * @param px Array of x positions.
 * @param py Array of y positions.
 * @param pz Array of z positions.
 * @param vx Array of x velocities.
 * @param vy Array of y velocities.
 * @param vz Array of z velocities.
 * @param dt Time step for the simulation.
 */
void SequentialAPAVXUpdate(const int n, double *m, double *px, double *py,
                           double *pz, double *vx, double *vy, double *vz,
                           const double dt) {
  static const __m256d DT = _mm256_set1_pd(dt);
  static const __m256d s = _mm256_set1_pd(_SOFTENING);

  for (int i = 0; i < n; ++i) {
    __m256d Fx = _mm256_set1_pd(0.0);
    __m256d Fy = _mm256_set1_pd(0.0);
    __m256d Fz = _mm256_set1_pd(0.0);

    for (int j = 0; j < n; j += 4) {  // exploit SIMD computing blocks of 4 pairs each time
      const __m256d Xi = _mm256_broadcast_sd(px + i);
      const __m256d Yi = _mm256_broadcast_sd(py + i);
      const __m256d Zi = _mm256_broadcast_sd(pz + i);

      const __m256d Xj = _mm256_loadu_pd(px + j);
      const __m256d Yj = _mm256_loadu_pd(py + j);
      const __m256d Zj = _mm256_loadu_pd(pz + j);

      const __m256d X = _mm256_sub_pd(Xj, Xi);
      const __m256d Y = _mm256_sub_pd(Yj, Yi);
      const __m256d Z = _mm256_sub_pd(Zj, Zi);

      __m256d M = _mm256_loadu_pd(m + j);

      const __m256d D = _mm256_fmadd_pd(X, X, _mm256_fmadd_pd(Y, Y, _mm256_fmadd_pd(Z, Z, _mm256_mul_pd(s, s))));
      const __m256d D_inv = fastinv(_mm256_sqrt_pd(D));
      const __m256d D_inv3 = _mm256_mul_pd(D_inv, _mm256_mul_pd(D_inv, D_inv));

      Fx = _mm256_fmadd_pd(D_inv3, _mm256_mul_pd(M, X), Fx);
      Fy = _mm256_fmadd_pd(D_inv3, _mm256_mul_pd(M, Y), Fy);
      Fz = _mm256_fmadd_pd(D_inv3, _mm256_mul_pd(M, Z), Fz);
    }
    const __m256d Vx = _mm256_mul_pd(Fx, DT);
    const __m256d Vy = _mm256_mul_pd(Fy, DT);
    const __m256d Vz = _mm256_mul_pd(Fz, DT);
    vx[i] += hsum_avx(Vx);
    vy[i] += hsum_avx(Vy);
    vz[i] += hsum_avx(Vz);
  }

  for (int i = 0; i < n; i += 4) {
    const __m256d X = _mm256_loadu_pd(px + i);
    const __m256d Y = _mm256_loadu_pd(py + i);
    const __m256d Z = _mm256_loadu_pd(pz + i);

    const __m256d Vx = _mm256_loadu_pd(vx + i);
    const __m256d Vy = _mm256_loadu_pd(vy + i);
    const __m256d Vz = _mm256_loadu_pd(vz + i);

    const __m256d Xn = _mm256_fmadd_pd(Vx, DT, X);
    const __m256d Yn = _mm256_fmadd_pd(Vy, DT, Y);
    const __m256d Zn = _mm256_fmadd_pd(Vz, DT, Z);

    _mm256_store_pd(px + i, Xn);
    _mm256_store_pd(py + i, Yn);
    _mm256_store_pd(pz + i, Zn);
  }
}
#endif

/**
 * @brief Simulate the n-body problem using AVX-optimized computations.
 * 
 * @tparam T Floating-point type (float or double).
 * @param n Number of bodies.
 * @param dt Time step for the simulation.
 * @param tEnd End time for the simulation.
 * @param seed Random seed for initializing bodies.
 */
template <typename T>
void SequentialAPAVXSimulate(int n, T dt, T tEnd, int seed) {
  T *m = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *px = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *py = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *pz = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *vx = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *vy = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *vz = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));

  // Initialize Bodies
  InitSoa<T>(n, m, px, py, pz, vx, vy, vz, seed);

#ifdef MONITOR_ENERGY
  T ek = EkSoa<T>(n, m, vx, vy, vz);
  T ep = EpSoa<T>(n, m, px, py, pz);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif

  TIMERSTART(simulation)
  // Simulation Loop
  for (T t = 0.0f; t < tEnd; t += dt) {
    // Update Bodies
    SequentialAPAVXUpdate(n, m, px, py, pz, vx, vy, vz, dt);
#ifdef MONITOR_ENERGY
    ek = EkSoa<T>(n, m, vx, vy, vz);
    ep = EpSoa<T>(n, m, px, py, pz);
    std::cout << "Etot: " << ek + ep << std::endl;
#endif
  }
#ifdef MONITOR_ENERGY
  ek = EkSoa<T>(n, m, vx, vy, vz);
  ep = EpSoa<T>(n, m, px, py, pz);
  std::cout << "Etot: " << ek + ep << std::endl;
#endif
  TIMERSTOP(simulation)

  // Free allocated memory
  _mm_free(m);
  _mm_free(px);
  _mm_free(py);
  _mm_free(pz);
  _mm_free(vx);
  _mm_free(vy);
  _mm_free(vz);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  int nbody = atoi(argv[1]);
  int seed = 0;
  if (argc == 3)
    seed = atoi(argv[2]);
#ifndef OMP
  srand(seed);
#endif
  SequentialAPAVXSimulate<MY_T>(nbody, 0.01, 1, seed);
}
