
#include "sequential_ap_avx.h"
#include <immintrin.h>// AVX intrinsics
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

// helper function given by SO - computes the horizontal sum
float hsum_avx(__m256 x) {
  // hiQuad = ( x7, x6, x5, x4 )
  const __m128 hiQuad = _mm256_extractf128_ps(x, 1);
  // loQuad = ( x3, x2, x1, x0 )
  const __m128 loQuad = _mm256_castps256_ps128(x);
  // sumQuad = ( x3 + x7, x2 + x6, x1 + x5, x0 + x4 )
  const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);
  // loDual = ( -, -, x1 + x5, x0 + x4 )
  const __m128 loDual = sumQuad;
  // hiDual = ( -, -, x3 + x7, x2 + x6 )
  const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);
  // sumDual = ( -, -, x1 + x3 + x5 + x7, x0 + x2 + x4 + x6 )
  const __m128 sumDual = _mm_add_ps(loDual, hiDual);
  // lo = ( -, -, -, x0 + x2 + x4 + x6 )
  const __m128 lo = sumDual;
  // hi = ( -, -, -, x1 + x3 + x5 + x7 )
  const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);
  // sum = ( -, -, -, x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 )
  const __m128 sum = _mm_add_ss(lo, hi);
  return _mm_cvtss_f32(sum);
}

void SequentialAPAVX::Update(float dt) {
  __m256 DT = _mm256_set_ps(dt, dt, dt, dt, dt, dt, dt, dt);
  __m256 s = _mm256_set_ps(1e-09f, 1e-09f, 1e-09f, 1e-09f, 1e-09f, 1e-09f, 1e-09f, 1e-09f);

  for (uint64_t i = 0; i < n_bodies; ++i) {

    __m256 Fx = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    __m256 Fy = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    __m256 Fz = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);

    for (uint64_t j = 0; j < n_bodies; j+=8) { // exploit the SIMD computing blocks of 8 pairs each time

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
      const __m256 mask = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f);
      M = _mm256_mul_ps(M, mask);

      const __m256 D = _mm256_fmadd_ps(X, X,
                                       _mm256_fmadd_ps(Y, Y,
                                                       _mm256_fmadd_ps(Z, Z, s))) ;
      const __m256 D_inv = _mm256_rsqrt_ps(D);
      const __m256 D_inv3 = _mm256_mul_ps(D_inv,
                                          _mm256_mul_ps(D_inv, D_inv)); //da scambiare con la riga precedente

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

  for (uint64_t i = 0; i < n_bodies; i+=8) {
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

std::ostream &operator<<(std::ostream &os, const SequentialAPAVX &nbody) {
  os << "Gravitational constant: " << nbody.G << std::endl;

  for (uint64_t i = 0; i < nbody.n_bodies; ++i) {
    os << "Body " << i + 1 << ":\n";
    os << "  Mass: " << nbody.m[i] << std::endl;
    os << "  Position (x, y, z): " << nbody.px[i] << ", " << nbody.py[i] << ", " << nbody.pz[i] << std::endl;
    os << "  Velocity (vx, vy, vz): " << nbody.vx[i] << ", " << nbody.vy[i] << ", " << nbody.vz[i] << std::endl;
  }
  return os;
}


void SequentialAPAVX::LogsToCSV(const std::string &filename) const {
    std::ofstream file(filename, std::ios_base::app);
    if (file.is_open()) {
      for (uint64_t i = 0; i < n_bodies; ++i) {
        file << px[i] << "," << py[i] << "," << pz[i];
        if (i < n_bodies - 1)
          file << ",";
      }
      file << std::endl;
      file.close();
    } else {
      std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
