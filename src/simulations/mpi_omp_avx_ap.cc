#include <cstring>
#include <iostream>
#include <omp.h>
#include <mpi.h>
#include "utilities/avx.h"
#include "utilities/time_utils.h"
#include "utilities/nbody_helpers.h"
#include "utilities/integrators.h"

// copyright NVIDIA
void SequentialAPAVXUpdate(const int n, const int localN, float *m, float *px, float *px_, float *py, float *py_, float *pz, float *pz_, float *vx, float *vx_, float *vy, float *vy_, float *vz, float *vz_, const float dt) {
  static __m256 DT = _mm256_set_ps(dt, dt, dt, dt, dt, dt, dt, dt);
  static __m256 s = _mm256_set_ps(
      _SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING,_SOFTENING);

  #pragma omp parallel for schedule(runtime)
  for (int i = 0; i < localN; ++i) {
    __m256 Fx = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    __m256 Fy = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    __m256 Fz = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);

    for (int j = 0; j < n; j += 8) {  // exploit the SIMD computing blocks of 8 pairs each time

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

  for (int i = 0; i < n; i += 8) {
    const __m256 X = _mm256_loadu_ps(px_ + i);
    const __m256 Y = _mm256_loadu_ps(py_ + i);
    const __m256 Z = _mm256_loadu_ps(pz_ + i);

    const __m256 Vx = _mm256_loadu_ps(vx_ + i);
    const __m256 Vy = _mm256_loadu_ps(vy_ + i);
    const __m256 Vz = _mm256_loadu_ps(vz_ + i);

    const __m256 Xn = _mm256_fmadd_ps(Vx, DT, X);
    const __m256 Yn = _mm256_fmadd_ps(Vy, DT, Y);
    const __m256 Zn = _mm256_fmadd_ps(Vz, DT, Z);

    _mm256_store_ps(px + i, Xn);
    _mm256_store_ps(pz + i, Zn);
    _mm256_store_ps(py + i, Yn);
  }
}


void MPIAPSimulate(uint64_t n, float dt, float tEnd, uint64_t seed){
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // we assume that n is a multiple of nproc
  int localN = n / nproc;

  float *m = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *px = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *py = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *pz = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *vx = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *vy = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));
  float *vz = static_cast<float *>(_mm_malloc(n * sizeof(float), 32));

  float *px_ = static_cast<float *>(_mm_malloc(localN * sizeof(float), 32));
  float *py_ = static_cast<float *>(_mm_malloc(localN * sizeof(float), 32));
  float *pz_ = static_cast<float *>(_mm_malloc(localN * sizeof(float), 32));
  float *vx_ = static_cast<float *>(_mm_malloc(localN * sizeof(float), 32));
  float *vy_ = static_cast<float *>(_mm_malloc(localN * sizeof(float), 32));
  float *vz_ = static_cast<float *>(_mm_malloc(localN * sizeof(float), 32));

  if (my_rank == 0)
    // Init Bodies
    InitSoa(n, m, px, py, pz, vx, vy, vz);

  MPI_Bcast(m, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(px, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(py, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(pz, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(vx, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(vy, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(vz, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

  for (int i = 0; i < localN; ++i) {
    int idx = my_rank*localN+i;
    px_[i] = px[idx];
    py_[i] = py[idx];
    pz_[i] = pz[idx];

    vx_[i] = vx[idx];
    vy_[i] = vy[idx];
    vz_[i] = vz[idx];
  }

  // Simulation Loop
  for (float t = 0.0f; t < tEnd; t += dt){
    // Update Bodies
    SequentialAPAVXUpdate(n, localN, m, px, px_, py, py_, pz, pz_, vx, vx_, vy, vy_, vz, vz_, dt);
    MPI_Allgather(px_, localN, MPI_FLOAT, px, localN, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(py_, localN, MPI_FLOAT, py, localN, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(pz_, localN, MPI_FLOAT, pz, localN, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(vx_, localN, MPI_FLOAT, vx, localN, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(vy_, localN, MPI_FLOAT, vy, localN, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(vz_, localN, MPI_FLOAT, vz, localN, MPI_FLOAT, MPI_COMM_WORLD);
  }
  if (my_rank == 0) {
    float ek = EkSoa(n, m, vx, vy, vz);
    float ep = EpSoa(n, m, px, py, pz);
    std::cout << "Etot: " << ek + ep << std::endl;
  }
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
  MPIAPSimulate(std::stoul(argv[1]), 0.01, 0.1, 0);
  TIMERSTOP(simulation)
  MPI_Finalize();
  //mpi_ap(12, 1.0f);
}