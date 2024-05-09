#include <iostream>
#include "simulations/sequential_ap_avx.h"
#include "utilities/avx.h"
#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

// copyright NVIDIA
void SequentialAPAVXUpdate(const uint64_t n, float *m, float *px, float *py, float *pz, float *vx, float *vy, float *vz, const float dt) {
  static __m256 DT = _mm256_set_ps(dt, dt, dt, dt, dt, dt, dt, dt);
  static __m256 s = _mm256_set_ps(
      _SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING);

  for (uint64_t i = 0; i < n; ++i) {
    __m256 Fx = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    __m256 Fy = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    __m256 Fz = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);

    for (uint64_t j = 0; j < n; j += 8) {  // exploit the SIMD computing blocks of 8 pairs each time

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
      //const __m256 mask =
      //   _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f);
      //M = _mm256_mul_ps(M, mask);

      const __m256 D = _mm256_fmadd_ps(
          X, X, _mm256_fmadd_ps(Y, Y, _mm256_fmadd_ps(Z, Z, _mm256_mul_ps(s,s))));
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


void SequentialAPAVXSimulate(uint64_t n, float dt, float tEnd, uint64_t seed){
  float *m = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *px = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *py = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *pz = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *vx = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *vy = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *vz = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));

  // Init Bodies
  InitSoa(n, m, px, py, pz, vx, vy, vz);

  // Simulation Loop
  for (float t = 0.0f; t < tEnd; t += dt){
    // Update Bodies
    SequentialAPAVXUpdate(n, m, px, py, pz, vx, vy, vz, dt);
    float ek = EkSoa(n, m, vx, vy ,vz);
    float ep = EpSoa(n, m, px, py, pz);
    std::cout << "Etot: " <<ek+ep <<std::endl;
  }
}

int main (int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  TIMERSTART(simulation)
  srand(0);
  SequentialAPAVXSimulate(std::stoul(argv[1]), 0.01, 0.1, 0);
  TIMERSTOP(simulation)
}