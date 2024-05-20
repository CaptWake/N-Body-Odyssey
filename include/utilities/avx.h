#ifndef UTILITIES_AVX_H_
#define UTILITIES_AVX_H_

#include <immintrin.h>

// helper function given by SO - computes the horizontal sum
#ifdef FLOAT
inline float hsum_avx(__m256 x) {
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
#else

inline double hsum_avx(__m256d x) {
  // hiBi = ( x3, x2 )
  const __m128d hiBi = _mm256_extractf128_pd(x, 1);
  // loBi = ( x1, x0 )
  const __m128d loBi = _mm256_castpd256_pd128(x);
  // sumBi = ( x3 + x1, x2 + x0 )
  const __m128d sumBi = _mm_add_pd(loBi, hiBi);
  // lo = ( -, x2 + x0 )
  const __m128d lo = sumBi;
  // hi = ( -, x3 + x1 )
  const __m128d hi = _mm_shuffle_pd(sumBi, sumBi, 0x1);
  // sum = ( -, x3 + x1 + x2 + x0 )
  const __m128d sum = _mm_add_sd(lo, hi);
  return _mm_cvtsd_f64(sum);
}
#endif

__m256d fastinv(__m256d y) {
  // exact results for powers of two
  __m256i const magic = _mm256_set1_epi64x(0x7fe0'0000'0000'0000);
  // Bit-magic: For powers of two this just inverts the exponent,
  // and values between that are linearly interpolated
  __m256d x =
      _mm256_castsi256_pd(_mm256_sub_epi64(magic, _mm256_castpd_si256(y)));

  // Newton-Raphson refinement: x = x*(2.0 - x*y):
  x = _mm256_mul_pd(x, _mm256_fnmadd_pd(x, y, _mm256_set1_pd(2.0)));

  return x;
}

#endif
