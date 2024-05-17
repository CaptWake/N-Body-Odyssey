#include "simulations/mpi_ap.h"

#include <mpi.h>

#include <cstring>
#include <iostream>

#include "utilities/integrators.h"
#include "utilities/nbody_helpers.h"
#include "utilities/time_utils.h"

void MPIAPUpdate(int localN, int n, const float *__restrict__ m,
                 const float *__restrict__ p, float *a,
                 float *__restrict__ m_rec, float *__restrict__ p_rec) {
  float *__restrict__ ai;
  float px, py, pz;
  long i, j;
  float dx, dy, dz, D;

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
      D = 1.0f / (D * sqrtf(D));
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
      D = 1.0f / (D * sqrtf(D));
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
  memcpy(m_rec, m, sizeof(float) * localN);
  memcpy(p_rec, p, sizeof(float) * localN * 3);

  for (int rounds = 1; rounds < world_size; ++rounds) {
    // fprintf(stderr, "Current rank: %d, src rank: %d, dest rank: %d\n",
    // world_rank, dst_rank, src_rank);
    MPI_Sendrecv_replace(m_rec,
                         localN,
                         MPI_FLOAT,
                         dst_rank,
                         0,
                         src_rank,
                         0,
                         MPI_COMM_WORLD,
                         &stat);
    MPI_Sendrecv_replace(p_rec,
                         localN * 3,
                         MPI_FLOAT,
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

static inline void performNBodyStep(const int localN, float *m, float *p,
                                    float *v, MPI_Request *requests,
                                    const float dt) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  int index;
  MPI_Status status;
  for (int n = 0; n < nproc; ++n) {
    MPI_Waitany(nproc, requests, &index, &status);
    int k = status.MPI_SOURCE;
    for (int i = my_rank * localN * 3; i < (my_rank + 1) * localN * 3; i += 3) {
      float fx = 0.0f;
      float fy = 0.0f;
      float fz = 0.0f;
      for (int j = k * localN * 3; j < (k + 1) * localN * 3; j += 3) {
        auto m2_id = j / 3;
        // compute distance pair
        auto dx = p[j] - p[i];
        auto dy = p[j + 1] - p[i + 1];
        auto dz = p[j + 2] - p[i + 2];

        auto d = dx * dx + dy * dy + dz * dz + _SOFTENING * _SOFTENING;
        auto d_inv = 1.0f / sqrtf(d);
        auto d_inv3 = d_inv * d_inv * d_inv;

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

void MPIAPSimulate(uint64_t n, float dt, float tEnd) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // we assume that n is a multiple of nproc
  int localN = n / nproc;
  float *pp = nullptr, *vv = nullptr, *aa = nullptr, *mm = nullptr;

  float *p = new float[localN * 3];
  float *v = new float[localN * 3];
  float *m = new float[localN];
  float *a = new float[localN * 3];
  float *p_rec = new float[localN * 3];
  float *m_rec = new float[localN * 3];

  if (my_rank == 0) {
    // Init Bodies
    mm = new float[n];
    pp = new float[3 * n];
    vv = new float[3 * n];
    aa = new float[3 * n];
    InitAos(n, mm, pp, vv, aa);
  }

  MPI_Scatter(
      pp, localN * 3, MPI_FLOAT, p, localN * 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(
      vv, localN * 3, MPI_FLOAT, v, localN * 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(
      aa, localN * 3, MPI_FLOAT, a, localN * 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(mm, localN, MPI_FLOAT, m, localN, MPI_FLOAT, 0, MPI_COMM_WORLD);

  // Simulation Loop
  for (float t = 0.0f; t < tEnd; t += dt) {
    performNBodyHalfStepA(localN, dt, p, v, a, m);
    // Update Bodies
    MPIAPUpdate(localN, n, m, p, a, m_rec, p_rec);
    performNBodyHalfStepB(localN, dt, p, v, a, m);
  }

  MPI_Gather(
      p, localN * 3, MPI_FLOAT, pp, localN * 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(
      v, localN * 3, MPI_FLOAT, vv, localN * 3, MPI_FLOAT, 0, MPI_COMM_WORLD);

  if (my_rank == 0) {
    float Epot = Ep(n, mm, pp);
    float Ekin = Ek(n, mm, vv);
    float E0 = Epot + Ekin;

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
void MPIAPSimulateV2(int n, float dt, float tEnd) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // we assume that n is a multiple of nproc
  int localN = n / nproc;
  float *m = new float[n];
  float *p = new float[3 * n];
  float *v = new float[3 * n];

  if (my_rank == 0)
    // Init Bodies
    InitAos(n, m, p, v);

  MPI_Request *requests = (MPI_Request *)malloc(nproc * sizeof(MPI_Request));

  MPI_Bcast(m, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(p,
              localN * 3,
              MPI_FLOAT,
              p + my_rank * localN * 3,
              localN * 3,
              MPI_FLOAT,
              0,
              MPI_COMM_WORLD);
  MPI_Scatter(v,
              localN * 3,
              MPI_FLOAT,
              v + my_rank * localN * 3,
              localN * 3,
              MPI_FLOAT,
              0,
              MPI_COMM_WORLD);

  int it = 0;
  // Simulation Loop
  for (float t = 0.0f; t < tEnd; t += dt) {
    // Update Bodies
    for (int i = 0; i < nproc; ++i) {
      MPI_Isend(p + my_rank * localN * 3,
                localN * 3,
                MPI_FLOAT,
                i,
                it,
                MPI_COMM_WORLD,
                &requests[i]);
      MPI_Irecv(p + i * localN * 3,
                localN * 3,
                MPI_FLOAT,
                i,
                it,
                MPI_COMM_WORLD,
                &requests[i]);
    }
    performNBodyStep(localN, m, p, v, requests, dt);
    ++it;
  }

  MPI_Gather(v + my_rank * localN * 3,
             localN * 3,
             MPI_FLOAT,
             v,
             3 * localN,
             MPI_FLOAT,
             0,
             MPI_COMM_WORLD);

  if (my_rank == 0) {
    float Epot = Ep(n, m, p);
    float Ekin = Ek(n, m, v);
    float E0 = Epot + Ekin;

    fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
    fprintf(stderr, "Eend: %.15g\n", E0);
  }

  delete[] m;
  delete[] p;
  delete[] v;
}

// Euler step https://en.wikipedia.org/wiki/File:Euler_leapfrog_comparison.gif//
void MPIAPSimulateV3(int n, float dt, float tEnd) {
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // we assume that n is a multiple of nproc
  int localN = n / nproc;
  float *m = new float[n];
  float *p = new float[3 * n];
  float *v = new float[3 * n];
  float *a = new float[3 * n];

  float *p_ = new float[localN * 3];
  float *v_ = new float[localN * 3];

  if (my_rank == 0)
    // Init Bodies
    InitAos(n, m, p, v, a);

  MPI_Bcast(m, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(p, 3 * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(v, 3 * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(a, 3 * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

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
    performNBodyStep(localN, n, m, p, p_, v, v_, dt);
    MPI_Allgather(
        p_, 3 * localN, MPI_FLOAT, p, 3 * localN, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(
        v_, 3 * localN, MPI_FLOAT, v, 3 * localN, MPI_FLOAT, MPI_COMM_WORLD);
  }

  if (my_rank == 0) {
    float Epot = Ep(n, m, p);
    float Ekin = Ek(n, m, v);
    float E0 = Epot + Ekin;

    fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
    fprintf(stderr, "Eend: %.15g\n", E0);
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
  int seed = 0;
  srand(seed);
  TIMERSTART(simulation)
  MPIAPSimulateV2(atoi(argv[1]), 0.01, 1);
  MPI_Barrier(MPI_COMM_WORLD);
  TIMERSTOP(simulation)
  MPI_Finalize();
}