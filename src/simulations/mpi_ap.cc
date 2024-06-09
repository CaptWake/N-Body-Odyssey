#include "simulations/mpi_ap.h"

#include <mpi.h>

#include <cstring>
#include <iostream>

#include "utilities/integrators.h"
#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"
#include <omp.h>

#ifdef DOUBLE
#define MPI_TYPE MPI_DOUBLE
#define MY_T double
#else
#define MPI_TYPE MPI_FLOAT
#define MY_T float
#endif

TIMERINIT(scatter)
TIMERINIT(broadcast)
TIMERINIT(gather)
TIMERINIT(allgather)
TIMERINIT(simulation)
TIMERINIT(waitany)
TIMERINIT(init)

template <typename T>
void MPIAPUpdate(int localN, int n, const T *__restrict__ m,
                 const T *__restrict__ p, T *a, T *__restrict__ m_rec,
                 T *__restrict__ p_rec) {
  T *__restrict__ ai;
  T px, py, pz;
  long i, j;
  T dx, dy, dz, D;

  // Force Calculation with Locals.
  for (i = 0; i < localN; ++i) {
    px = p[3 * i + 0];
    py = p[3 * i + 1];
    pz = p[3 * i + 2];
    ai = a + 3 * i;
    ai[0] = 0.0;
    ai[1] = 0.0;
    ai[2] = 0.0;

    for (j = 0; j < i; ++j) {
      // really dx is other way around, but this way we can avoid -1.0* later.
      dx = p[3 * j] - px;
      dy = p[3 * j + 1] - py;
      dz = p[3 * j + 2] - pz;
      D = dx * dx + dy * dy + dz * dz;
      D += _SOFTENING * _SOFTENING;
      D = 1.0 / (D * sqrt(D));
      ai[0] += m[j] * dx * D;
      ai[1] += m[j] * dy * D;
      ai[2] += m[j] * dz * D;
    }
    for (j = i + 1; j < localN; ++j) {
      dx = p[3 * j] - px;
      dy = p[3 * j + 1] - py;
      dz = p[3 * j + 2] - pz;
      D = dx * dx + dy * dy + dz * dz;
      D += _SOFTENING * _SOFTENING;
      D = 1.0 / (D * sqrt(D));
      ai[0] += m[j] * dx * D;
      ai[1] += m[j] * dy * D;
      ai[2] += m[j] * dz * D;
    }
  }
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int src_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;
  int dst_rank = (world_rank + 1) % world_size;
  MPI_Status stat;

  // fill send buffers initially
  memcpy(m_rec, m, sizeof(T) * localN);
  memcpy(p_rec, p, sizeof(T) * localN * 3);

  for (int rounds = 1; rounds < world_size; ++rounds) {
    // fprintf(stderr, "Current rank: %d, src rank: %d, dest rank: %d\n",
    // world_rank, dst_rank, src_rank);
    MPI_Sendrecv_replace(m_rec,
                         localN,
                         MPI_TYPE,
                         dst_rank,
                         0,
                         src_rank,
                         0,
                         MPI_COMM_WORLD,
                         &stat);
    MPI_Sendrecv_replace(p_rec,
                         localN * 3,
                         MPI_TYPE,
                         dst_rank,
                         0,
                         src_rank,
                         0,
                         MPI_COMM_WORLD,
                         &stat);

    // update each local body with the influence from the new bodies
    //  Force Calculation with Remotes
    // TODO could use blocking here for locality
    for (i = 0; i < localN; ++i) {
      px = p[3 * i + 0];
      py = p[3 * i + 1];
      pz = p[3 * i + 2];
      ai = a + 3 * i;

      for (j = 0; j < localN; ++j) {
        // really dx is other way around, but this way we can avoid -1.0* later.
        dx = p_rec[3 * j] - px;
        dy = p_rec[3 * j + 1] - py;
        dz = p_rec[3 * j + 2] - pz;
        D = dx * dx + dy * dy + dz * dz;
        D += _SOFTENING * _SOFTENING;
        D = 1.0 / (D * sqrt(D));
        ai[0] += m_rec[j] * dx * D;
        ai[1] += m_rec[j] * dy * D;
        ai[2] += m_rec[j] * dz * D;
      }
    }
  }
}

template <typename T>
static inline void performNBodyStep(const int localN, T *m, T *p, T *v,
                                    MPI_Request *requests, const T dt) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  int index;
  MPI_Status status;
  for (int n = 0; n < nproc; ++n) {
    TIMERSTART(waitany)
    MPI_Waitany(nproc, requests, &index, &status);
    int k = status.MPI_SOURCE;
    for (int i = my_rank * localN * 3; i < (my_rank + 1) * localN * 3; i += 3) {
      T fx = 0.0;
      T fy = 0.0;
      T fz = 0.0;
      for (int j = k * localN * 3; j < (k + 1) * localN * 3; j += 3) {
        uint64_t m2_id = j / 3;
        // compute distance pair
        T dx = p[j] - p[i];
        T dy = p[j + 1] - p[i + 1];
        T dz = p[j + 2] - p[i + 2];

        T d = dx * dx + dy * dy + dz * dz + _SOFTENING * _SOFTENING;
        T d_inv = 1.0 / sqrt(d);
        T d_inv3 = d_inv * d_inv * d_inv;

        fx += d_inv3 * m[m2_id] * dx;
        fy += d_inv3 * m[m2_id] * dy;
        fz += d_inv3 * m[m2_id] * dz;
      }

      v[i] += fx * dt;
      v[i + 1] += fy * dt;
      v[i + 2] += fz * dt;
    }
  }
  for (int i = my_rank * localN * 3; i < (my_rank + 1) * localN * 3; i += 3) {
    p[i] += v[i] * dt;
    p[i + 1] += v[i + 1] * dt;
    p[i + 2] += v[i + 2] * dt;
  }
}

template <typename T>
void MPIAPSimulate(uint64_t n, T dt, T tEnd) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // we assume that n is a multiple of nproc
  int localN = n / nproc;
  T *pp = nullptr, *vv = nullptr, *aa = nullptr, *mm = nullptr;

  T *p = new T[localN * 3];
  T *v = new T[localN * 3];
  T *m = new T[localN];
  T *a = new T[localN * 3];
  T *p_rec = new T[localN * 3];
  T *m_rec = new T[localN * 3];

  if (my_rank == 0) {
    // Init Bodies
    mm = new T[n];
    pp = new T[3 * n];
    vv = new T[3 * n];
    aa = new T[3 * n];
    InitAos<T>(n, mm, pp, vv, aa);
  }

  TIMERSTART(scatter)
  MPI_Scatter(
      pp, localN * 3, MPI_TYPE, p, localN * 3, MPI_TYPE, 0, MPI_COMM_WORLD);
  MPI_Scatter(
      vv, localN * 3, MPI_TYPE, v, localN * 3, MPI_TYPE, 0, MPI_COMM_WORLD);
  MPI_Scatter(
      aa, localN * 3, MPI_TYPE, a, localN * 3, MPI_TYPE, 0, MPI_COMM_WORLD);
  MPI_Scatter(mm, localN, MPI_TYPE, m, localN, MPI_TYPE, 0, MPI_COMM_WORLD);
  TIMERSTOP(scatter)

  // Simulation Loop
  TIMERSTART(simulation)
  for (T t = 0.0f; t < tEnd; t += dt) {
    performNBodyHalfStepA<T>(localN, dt, p, v, a, m);
    // Update Bodies
    MPIAPUpdate<T>(localN, n, m, p, a, m_rec, p_rec);
    performNBodyHalfStepB<T>(localN, dt, p, v, a, m);
  }
  TIMERSTOP(simulation)

  TIMERSTART(gather)
  MPI_Gather(
      p, localN * 3, MPI_TYPE, pp, localN * 3, MPI_TYPE, 0, MPI_COMM_WORLD);
  MPI_Gather(
      v, localN * 3, MPI_TYPE, vv, localN * 3, MPI_TYPE, 0, MPI_COMM_WORLD);
  TIMERSTOP(gather)

  TIMERPRINT(scatter)
  TIMERPRINT(gather)
  TIMERPRINT(simulation)

  MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == 0) {
    T Epot = Ep<T>(n, mm, pp);
    T Ekin = Ek<T>(n, mm, vv);
    T E0 = Epot + Ekin;

    fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
    fprintf(stderr, "Eend: %.15g\n", E0);

    delete[] pp;
    delete[] vv;
    delete[] aa;
    delete[] mm;
  }

  delete[] p;
  delete[] v;
  delete[] a;
  delete[] m;
  delete[] p_rec;
  delete[] m_rec;
}

// Euler step https://en.wikipedia.org/wiki/File:Euler_leapfrog_comparison.gif//
template <typename T>
void MPIAPSimulateV2(int n, T dt, T tEnd, int seed) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // we assume that n is a multiple of nproc
  int localN = n / nproc;
  T *m = new T[n];
  T *p = new T[3 * n];
  T *v = new T[3 * n];

  if (my_rank == 0){
    // Init Bodies
    TIMERSTART(init)
    InitAos<T>(n, m, p, v, seed);
    TIMERSTOP(init)
    TIMERPRINT(init)
   }

  MPI_Request *requests = (MPI_Request *)malloc(nproc * sizeof(MPI_Request));

  TIMERSTART(simulation)
  TIMERSTART(broadcast)
  MPI_Bcast(m, n, MPI_TYPE, 0, MPI_COMM_WORLD);
  TIMERSTOP(broadcast)
  TIMERPRINT(broadcast)

  TIMERSTART(scatter)
  MPI_Scatter(p,
              localN * 3,
              MPI_TYPE,
              p + my_rank * localN * 3,
              localN * 3,
              MPI_TYPE,
              0,
              MPI_COMM_WORLD);
  MPI_Scatter(v,
              localN * 3,
              MPI_TYPE,
              v + my_rank * localN * 3,
              localN * 3,
              MPI_TYPE,
              0,
              MPI_COMM_WORLD);
  TIMERSTOP(scatter)
  TIMERPRINT(scatter)

  int it = 0;
  // Simulation Loop
  for (T t = 0.0f; t < tEnd; t += dt) {
    // Update Bodies
    for (int i = 0; i < nproc; ++i) {
      TIMERSTART(waitany)
      MPI_Isend(p + my_rank * localN * 3,
                localN * 3,
                MPI_TYPE,
                i,
                it,
                MPI_COMM_WORLD,
                &requests[i]);
      MPI_Irecv(p + i * localN * 3,
                localN * 3,
                MPI_TYPE,
                i,
                it,
                MPI_COMM_WORLD,
                &requests[i]);
       TIMERSTOP(waitany)
    }
    performNBodyStep<T>(localN, m, p, v, requests, dt);
    ++it;
  }
  TIMERSTOP(simulation)
  TIMERPRINT(waitany)
  TIMERPRINT(simulation)

  TIMERSTART(gather)
  MPI_Gather(v + my_rank * localN * 3,
             localN * 3,
             MPI_TYPE,
             v,
             3 * localN,
             MPI_TYPE,
             0,
             MPI_COMM_WORLD);
  TIMERSTOP(gather)
  TIMERPRINT(gather)

  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0) {
    T Epot = Ep<T>(n, m, p);
    T Ekin = Ek<T>(n, m, v);
    T E0 = Epot + Ekin;

    fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
    fprintf(stderr, "Eend: %.15g\n", E0);
  }

  delete[] m;
  delete[] p;
  delete[] v;
}

// Euler step https://en.wikipedia.org/wiki/File:Euler_leapfrog_comparison.gif//
template <typename T>
void MPIAPSimulateV3(int n, T dt, T tEnd, int seed) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // we assume that n is a multiple of nproc
  int localN = n / nproc;
  T *m = new T[n];
  T *p = new T[3 * n];
  T *v = new T[3 * n];
  T *a = new T[3 * n];

  T *p_ = new T[localN * 3];
  T *v_ = new T[localN * 3];

  if (my_rank == 0)
    // Init Bodies
    InitAos<T>(n, m, p, v, a, seed);

  TIMERSTART(simulation)
  TIMERSTART(broadcast)
  MPI_Bcast(m, n, MPI_TYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(p, 3 * n, MPI_TYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(v, 3 * n, MPI_TYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(a, 3 * n, MPI_TYPE, 0, MPI_COMM_WORLD);
  TIMERSTOP(broadcast)

  for (int i = 0; i < localN * 3; i += 3) {
    p_[i] = p[my_rank * localN * 3 + i];
    p_[i + 1] = p[my_rank * localN * 3 + i + 1];
    p_[i + 2] = p[my_rank * localN * 3 + i + 2];

    v_[i] = v[my_rank * localN * 3 + i];
    v_[i + 1] = v[my_rank * localN * 3 + i + 1];
    v_[i + 2] = v[my_rank * localN * 3 + i + 2];
  }

  // Simulation Loop
  for (float t = 0.0f; t < tEnd; t += dt) {
    // Update Bodies
    performNBodyStep<T>(localN, n, m, p, p_, v, v_, dt);
    TIMERSTART(allgather)
    MPI_Allgather(
        p_, 3 * localN, MPI_TYPE, p, 3 * localN, MPI_TYPE, MPI_COMM_WORLD);
    MPI_Allgather(
        v_, 3 * localN, MPI_TYPE, v, 3 * localN, MPI_TYPE, MPI_COMM_WORLD);
    TIMERSTOP(allgather)

  }
  TIMERSTOP(simulation)

  TIMERPRINT(broadcast)
  TIMERPRINT(simulation)
  TIMERPRINT(allgather)

  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0) {
    T Epot = Ep<T>(n, m, p);
    T Ekin = Ek<T>(n, m, v);
    T E0 = Epot + Ekin;

    printf("Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
    printf("Eend: %.15g\n", E0);
  }
  delete[] m;
  delete[] p;
  delete[] v;
  delete[] a;
  delete[] p_;
  delete[] v_;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
#ifndef OMP
  srand(atoi(argv[2]));
#endif
  MPIAPSimulateV2<MY_T>(atoi(argv[1]), 0.01, 10, atoi(argv[2]));
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
