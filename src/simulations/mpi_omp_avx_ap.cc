#include <mpi.h>
#include <omp.h>

#include <cstring>
#include <iostream>

#include "utilities/avx.h"
#include "utilities/integrators.h"
#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

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

void MPIAPSimulate(uint64_t n, float dt, float tEnd, uint64_t seed) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  // we assume that n is a multiple of nproc
  int localN = n / nproc;

  MPI_Datatype block;
  MPI_Type_vector(6, localN, n, MPI_FLOAT, &block);
  MPI_Type_commit(&block);

  float *m = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *pv = static_cast<float *>(_mm_malloc(6 * n * sizeof(float), 32));
  float *pv_ = static_cast<float *>(_mm_malloc(6 * localN * sizeof(float), 32));

  if (my_rank == 0)
    // Init Bodies
    InitSoa(n, m, pv, pv + n, pv + 2 * n, pv + 3 * n, pv + 4 * n, pv + 5 * n);

  MPI_Bcast(m, 3 * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(pv, 6 * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

  for (int i = 0; i < localN; ++i) {
    int idx = my_rank * localN + i;
    pv_[i] = pv[idx];
    pv_[i + localN] = pv[idx + n];
    pv_[i + 2 * localN] = pv[idx + 2 * n];
    pv_[i + 3 * localN] = pv[idx + 3 * n];
    pv_[i + 4 * localN] = pv[idx + 4 * n];
    pv_[i + 5 * localN] = pv[idx + 5 * n];
  }

  // Simulation Loop
  for (float t = 0.0f; t < tEnd; t += dt) {
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
        memcpy(&pv[i * n], &pv_[i * localN], localN * sizeof(float));

      for (int i = 1; i < nproc; i++)
        MPI_Recv(
            &pv[i * localN], 1, block, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Send(pv_, 6 * localN, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(pv, 6 * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }
  if (my_rank == 0) {
    float ek = EkSoa(n, m, pv + 3 * n, pv + 4 * n, pv + 5 * n);
    float ep = EpSoa(n, m, pv, pv + n, pv + 2 * n);
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
  MPIAPSimulate(std::stoul(argv[1]), 0.01, 1, 0);
  MPI_Barrier(MPI_COMM_WORLD);
  TIMERSTOP(simulation)
  MPI_Finalize();
  // mpi_ap(12, 1.0f);
}