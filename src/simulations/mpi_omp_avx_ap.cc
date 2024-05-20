#include <mpi.h>
#include <omp.h>

#include <cstring>
#include <iostream>

#include "utilities/avx.h"
#include "utilities/integrators.h"
#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

#ifdef DOUBLE
#define MPI_TYPE MPI_DOUBLE
#define MY_T double
#else
#define MPI_TYPE MPI_FLOAT
#define MY_T float
#endif


#ifdef FLOAT
// copyright NVIDIA
void SequentialAPAVXUpdate(const int n, const int localN, float *m, float *px,
                           float *px_, float *py, float *py_, float *pz,
                           float *pz_, float *vx, float *vx_, float *vy,
                           float *vy_, float *vz, float *vz_, const float dt) {
  static __m256 DT = _mm256_set1_ps(dt);
  static __m256 s = _mm256_set1_ps(_SOFTENING);

#pragma omp parallel for schedule(runtime)
  for (int i = 0; i < localN; ++i) {
    __m256 Fx = _mm256_set1_ps(0.0f);
    __m256 Fy = _mm256_set1_ps(0.0f);
    __m256 Fz = _mm256_set1_ps(0.0f);

    for (int j = 0; j < n;
         j += 8) {  // exploit the SIMD computing blocks of 8 pairs each time

      const __m256 Xi = _mm256_broadcast_ss(px_ + i);
      const __m256 Yi = _mm256_broadcast_ss(py_ + i);
      const __m256 Zi = _mm256_broadcast_ss(pz_ + i);

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
    vx_[i] += hsum_avx(Vx);
    vy_[i] += hsum_avx(Vy);
    vz_[i] += hsum_avx(Vz);
  }

  for (int i = 0; i < localN; i += 8) {
    const __m256 X = _mm256_loadu_ps(px_ + i);
    const __m256 Y = _mm256_loadu_ps(py_ + i);
    const __m256 Z = _mm256_loadu_ps(pz_ + i);

    const __m256 Vx = _mm256_loadu_ps(vx_ + i);
    const __m256 Vy = _mm256_loadu_ps(vy_ + i);
    const __m256 Vz = _mm256_loadu_ps(vz_ + i);

    const __m256 Xn = _mm256_fmadd_ps(Vx, DT, X);
    const __m256 Yn = _mm256_fmadd_ps(Vy, DT, Y);
    const __m256 Zn = _mm256_fmadd_ps(Vz, DT, Z);

    _mm256_store_ps(px_ + i, Xn);
    _mm256_store_ps(py_ + i, Yn);
    _mm256_store_ps(pz_ + i, Zn);
  }
}

#else

// copyright NVIDIA
void SequentialAPAVXUpdate(const int n, const int localN, double *m, double *px,
                           double *px_, double *py, double *py_, double *pz,
                           double *pz_, double *vx, double *vx_, double *vy,
                           double *vy_, double *vz, double *vz_, const double dt) {
  static __m256d DT = _mm256_set1_pd(dt);
  static __m256d s = _mm256_set1_pd(_SOFTENING);

#pragma omp parallel for schedule(runtime)
  for (int i = 0; i < localN; ++i) {
    __m256d Fx = _mm256_set1_pd(0.0);
    __m256d Fy = _mm256_set1_pd(0.0);
    __m256d Fz = _mm256_set1_pd(0.0);

    for (int j = 0; j < n;
         j += 4) {  // exploit the SIMD computing blocks of 8 pairs each time

      const __m256d Xi = _mm256_broadcast_sd(px_ + i);
      const __m256d Yi = _mm256_broadcast_sd(py_ + i);
      const __m256d Zi = _mm256_broadcast_sd(pz_ + i);

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
    vx_[i] += hsum_avx(Vx);
    vy_[i] += hsum_avx(Vy);
    vz_[i] += hsum_avx(Vz);
  }

  for (int i = 0; i < localN; i += 4) {
    const __m256d X = _mm256_loadu_pd(px_ + i);
    const __m256d Y = _mm256_loadu_pd(py_ + i);
    const __m256d Z = _mm256_loadu_pd(pz_ + i);

    const __m256d Vx = _mm256_loadu_pd(vx_ + i);
    const __m256d Vy = _mm256_loadu_pd(vy_ + i);
    const __m256d Vz = _mm256_loadu_pd(vz_ + i);

    const __m256d Xn = _mm256_fmadd_pd(Vx, DT, X);
    const __m256d Yn = _mm256_fmadd_pd(Vy, DT, Y);
    const __m256d Zn = _mm256_fmadd_pd(Vz, DT, Z);

    _mm256_store_pd(px_ + i, Xn);
    _mm256_store_pd(py_ + i, Yn);
    _mm256_store_pd(pz_ + i, Zn);
  }
}


#endif
template<typename T>
void MPIAPSimulate(uint64_t n, T dt, T tEnd, uint64_t seed) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  // we assume that n is a multiple of nproc
  int localN = n / nproc;

  MPI_Datatype block;
  MPI_Type_vector(6, localN, n, MPI_TYPE, &block);
  MPI_Type_commit(&block);

  T *m = static_cast<T *>(_mm_malloc(n * sizeof(T), 32));
  T *pv = static_cast<T *>(_mm_malloc(6 * n * sizeof(T), 32));
  T *pv_ = static_cast<T *>(_mm_malloc(6 * localN * sizeof(T), 32));


  if (my_rank == 0)
    // Init Bodies
    InitSoa<T>(n, m, pv, pv + n, pv + 2 * n, pv + 3 * n, pv + 4 * n, pv + 5 * n);

  MPI_Bcast(m, n, MPI_TYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(pv, 6 * n, MPI_TYPE, 0, MPI_COMM_WORLD);

  for (int i = 0; i < localN; ++i) {
    int idx = my_rank * localN + i;
    pv_[i] = pv[idx];
    pv_[i + localN] = pv[idx + n];
    pv_[i + 2 * localN] = pv[idx + 2 * n];
    pv_[i + 3 * localN] = pv[idx + 3 * n];
    pv_[i + 4 * localN] = pv[idx + 4 * n];
    pv_[i + 5 * localN] = pv[idx + 5 * n];
  }

  // Simulation Loopp
  for (T t = 0.0f; t < tEnd; t += dt) {
    // Update Bodies
    SequentialAPAVXUpdate(n,
                          localN,
                          m,
                          pv,
                          pv_,
                          pv + n,
                          pv_ + localN,
                          pv + 2 * n,
                          pv_ + 2 * localN,
                          pv + 3 * n,
                          pv_ + 3 * localN,
                          pv + 4 * n,
                          pv_ + 4 * localN,
                          pv + 5 * n,
                          pv_ + 5 * localN,
                          dt);

    if (my_rank == 0) {
      for (int i = 0; i < 6; i++)
        memcpy(&pv[i * n], &pv_[i * localN], localN * sizeof(T));

      for (int i = 1; i < nproc; i++)
        MPI_Recv(
            &pv[i * localN], 1, block, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Send(pv_, 6 * localN, MPI_TYPE, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(pv, 6 * n, MPI_TYPE, 0, MPI_COMM_WORLD);
  }
  if (my_rank == 0) {
    T ek = EkSoa<T>(n, m, pv + 3 * n, pv + 4 * n, pv + 5 * n);
    T ep = EpSoa<T>(n, m, pv, pv + n, pv + 2 * n);
    std::cout << "Etot: " << ek + ep << std::endl;
  }
  MPI_Type_free(&block);  // Free the datatype when done
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  srand(0);
  SetScheduleType(argv[2], atoi(argv[3]));
  SetNumThread(atoi(argv[4]));
  TIMERSTART(simulation)
  MPIAPSimulate<MY_T>(std::stoul(argv[1]), 0.01, 1, 0);
  MPI_Barrier(MPI_COMM_WORLD);
  TIMERSTOP(simulation)
  MPI_Finalize();
  // mpi_ap(12, 1.0f);
}