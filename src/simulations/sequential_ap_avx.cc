#include "simulations/sequential_ap_avx.h"

#include <iostream>

#include "utilities/avx.h"
#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

// copyright NVIDIA

#ifdef FLOAT
void SequentialAPAVXUpdate(const uint64_t n, float *m, float *px, float *py,
                           float *pz, float *vx, float *vy, float *vz,
                           const float dt) {
  static __m256 DT = _mm256_set1_ps(dt);
  static __m256 s = _mm256_set1_ps(_SOFTENING);

  for (uint64_t i = 0; i < n; ++i) {
    __m256 Fx = _mm256_set1_ps(0.0f);
    __m256 Fy = _mm256_set1_ps(0.0f);
    __m256 Fz = _mm256_set1_ps(0.0f);

    for (uint64_t j = 0; j < n;
         j += 8) {  // exploit the SIMD computing blocks of 8 pairs each time

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
      // const __m256 mask =
      //    _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f);
      // M = _mm256_mul_ps(M, mask);

      const __m256 D = _mm256_fmadd_ps(
          X,
          X,
          _mm256_fmadd_ps(Y, Y, _mm256_fmadd_ps(Z, Z, _mm256_mul_ps(s, s))));
      const __m256 D_inv = _mm256_rsqrt_ps(D);
      const __m256 D_inv3 = _mm256_mul_ps(
          D_inv,
          _mm256_mul_ps(D_inv, D_inv));  // da scambiare con la riga precedente

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

  for (uint64_t i = 0; i < n; i += 8) {
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
    _mm256_store_ps(pz + i, Zn);
    _mm256_store_ps(py + i, Yn);
  }
}
#else
void SequentialAPAVXUpdate(const uint64_t n, double *m, double *px, double *py,
                           double *pz, double *vx, double *vy, double *vz,
                           const double dt) {
  static const __m256d DT = _mm256_set1_pd(dt);
  static const __m256d s = _mm256_set1_pd(_SOFTENING);

  for (uint64_t i = 0; i < n; ++i) {
    __m256d Fx = _mm256_set1_pd(0.0);
    __m256d Fy = _mm256_set1_pd(0.0);
    __m256d Fz = _mm256_set1_pd(0.0);

    for (uint64_t j = 0; j < n;
         j += 4) {  // exploit the SIMD computing blocks of 8 pairs each time

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
      // const __m256 mask =
      //    _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f);
      // M = _mm256_mul_ps(M, mask);

      const __m256d D = _mm256_fmadd_pd(
          X,
          X,
          _mm256_fmadd_pd(Y, Y, _mm256_fmadd_pd(Z, Z, _mm256_mul_pd(s, s))));
      const __m256d D_sqrt = _mm256_sqrt_pd(D);
      const __m256d D_inv = fastinv(D_sqrt);
      const __m256d D_inv3 = _mm256_mul_pd(
          D_inv,
          _mm256_mul_pd(D_inv, D_inv));  // da scambiare con la riga precedente

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

  for (uint64_t i = 0; i < n; i += 4) {
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
    _mm256_store_pd(pz + i, Zn);
    _mm256_store_pd(py + i, Yn);
  }
}
#endif

template <typename T>
void SequentialAPAVXSimulate(uint64_t n, T dt, T tEnd, uint64_t seed) {
  T *m = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *px = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *py = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *pz = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *vx = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *vy = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *vz = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));

  // Init Bodies
  InitSoa<T>(n, m, px, py, pz, vx, vy, vz);

  TIMERSTART(simulation)
  // Simulation Loop
  for (T t = 0.0f; t < tEnd; t += dt) {
    // Update Bodies
    SequentialAPAVXUpdate(n, m, px, py, pz, vx, vy, vz, dt);
  }
  TIMERSTOP(simulation)

  float ek = EkSoa<T>(n, m, vx, vy, vz);
  float ep = EpSoa<T>(n, m, px, py, pz);
  std::cout << "Etot: " << ek + ep << std::endl;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  srand(0);
  SequentialAPAVXSimulate<MY_T>(std::stoul(argv[1]), 0.01, 1, 0);
}