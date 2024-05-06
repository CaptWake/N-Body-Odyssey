#include "mpi_ap.h"

#include <cstring>

#define SOFTENING 1e-9

void computeForcesNaive_NBMPI(int localN,
                              const float* __restrict__ m,
                              const float* __restrict__ p,
                              float* a,
                              float* __restrict__ m_rec,
                              float* __restrict__ p_rec)
{
  float* __restrict__ ai;
  float px, py, pz;
  long i, j;
  float dx, dy, dz, D;

  for (i = 0; i < localN; ++i) {
    px = p[3*i + 0];
    py = p[3*i + 1];
    pz = p[3*i + 2];
    ai = a + 3*i;
    ai[0] = 0.0; ai[1] = 0.0; ai[2] = 0.0;

    for (j = 0; j < i; ++j) {
      //really dx is other way around, but this way we can avoid -1.0* later.
      dx = p[3*j] - px;
      dy = p[3*j + 1] - py;
      dz = p[3*j + 2] - pz;
      D = dx*dx + dy*dy + dz*dz;
      D += SOFTENING*SOFTENING;
      D = 1.0f / (D*sqrtf(D));
      ai[0] += m[j]*dx*D;
      ai[1] += m[j]*dy*D;
      ai[2] += m[j]*dz*D;
    }
    for (j = i+1; j < localN; ++j) {
      dx = p[3*j] - px;
      dy = p[3*j + 1] - py;
      dz = p[3*j + 2] - pz;
      D = dx*dx + dy*dy + dz*dz;
      D += SOFTENING*SOFTENING;
      D = 1.0f / (D*sqrtf(D));
      ai[0] += m[j]*dx*D;
      ai[1] += m[j]*dy*D;
      ai[2] += m[j]*dz*D;
    }
  }
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int src_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;
  int dst_rank = (world_rank + 1) % world_size;
  MPI_Status stat;

  //fill send buffers initially
  memcpy(m_rec, m, sizeof(float)*localN);
  memcpy(p_rec, p, sizeof(float)*localN*3);

  for (int rounds = 1; rounds < world_size; ++rounds) {
    // fprintf(stderr, "Current rank: %d, src rank: %d, dest rank: %d\n", world_rank, dst_rank, src_rank);
    MPI_Sendrecv_replace(m_rec, localN, MPI_FLOAT,
                         dst_rank, 0,
                         src_rank, 0,
                         MPI_COMM_WORLD, &stat);
    MPI_Sendrecv_replace(p_rec, localN*3, MPI_FLOAT,
                         dst_rank, 0,
                         src_rank, 0,
                         MPI_COMM_WORLD, &stat);

    //update each local body with the influence from the new bodies
    //TODO could use blocking here for locality
    for (i = 0; i < localN; ++i) {
      px = p[3*i + 0];
      py = p[3*i + 1];
      pz = p[3*i + 2];
      ai = a + 3*i;

      for (j = 0; j < localN; ++j) {
        //really dx is other way around, but this way we can avoid -1.0* later.
        dx = p_rec[3*j] - px;
        dy = p_rec[3*j + 1] - py;
        dz = p_rec[3*j + 2] - pz;
        D = dx*dx + dy*dy + dz*dz;
        D += SOFTENING*SOFTENING;
        D = 1.0 / (D*sqrt(D));
        ai[0] += m_rec[j]*dx*D;
        ai[1] += m_rec[j]*dy*D;
        ai[2] += m_rec[j]*dz*D;
      }
    }
  }
}

void mpi_ap(int n_bodies, const float grav_const) {
  MPI_Init(nullptr, nullptr);
  static std::random_device rd;  // random device engine, usually based on
                                 // /dev/random on UNIX-like systems
  float dt = 0.01f;
  int my_rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // we assume that the number of bodies is a multiple of the number of processes
  int local_n_bodies = n_bodies / nproc;
  float *p = nullptr, *v = nullptr, *a = nullptr, *m = nullptr, *p_rec = nullptr, *m_rec = nullptr;
  float *pp = nullptr, *vv = nullptr, *aa = nullptr, *mm = nullptr;
  if (my_rank == 0) {
    // initialize Mersennes' twister using rd to generate the seed
    static std::mt19937 engine{0};  // rd()};
    std::uniform_real_distribution<float> density(-1, 1);
    const uint64_t n_coords = n_bodies * 3;

    pp = new float[n_bodies*3];
    vv = new float[n_bodies*3];
    aa = new float[n_bodies*3];
    mm = new float[n_bodies];

    for (uint64_t i = 0; i < n_bodies*3; i+=3) {
      pp[i] = i;//density(engine);
      pp[i+1] = i+1;//density(engine);
      pp[i+2] = i+2;//density(engine);

      vv[i] = i;//density(engine);
      vv[i+1] = i+1; //density(engine);
      vv[i+2] = i+2; //density(engine);

      mm[i/3] = i/3;//density(engine);
    }
  }

  p = new float[local_n_bodies * 3];
  v = new float[local_n_bodies * 3];
  m = new float[local_n_bodies * 3];
  a = new float[local_n_bodies * 3];
  p_rec = new float[local_n_bodies * 3];
  m_rec = new float[local_n_bodies * 3];

  MPI_Scatter(pp, local_n_bodies*3, MPI_FLOAT, p, local_n_bodies*3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(vv, local_n_bodies*3, MPI_FLOAT, v, local_n_bodies*3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(aa, local_n_bodies*3, MPI_FLOAT, a, local_n_bodies*3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(mm, local_n_bodies, MPI_FLOAT, m, local_n_bodies, MPI_FLOAT, 0, MPI_COMM_WORLD);

  for (long i = 0; i < local_n_bodies; ++i) {
    //kick, drift
    v[3*i + 0] += 0.5f * a[3*i + 0] * dt;
    p[3*i + 0] += v[3*i + 0] * dt;
    v[3*i + 1] += 0.5f * a[3*i + 1] * dt;
    p[3*i + 1] += v[3*i + 1] * dt;
    v[3*i + 2] += 0.5f * a[3*i + 2] * dt;
    p[3*i + 2] += v[3*i + 2] * dt;
  }

  // update acceleration
  computeForcesNaive_NBMPI(local_n_bodies, m, p, a, m_rec, p_rec);

  for (long i = 0; i < local_n_bodies; ++i) {
    //kick
    v[3*i + 0] += 0.5f * a[3*i + 0] * dt;
    v[3*i + 1] += 0.5f * a[3*i + 1] * dt;
    v[3*i + 2] += 0.5f * a[3*i + 2] * dt;
  }

  MPI_Gather(p, local_n_bodies*3, MPI_FLOAT, pp, local_n_bodies*3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(v, local_n_bodies*3, MPI_FLOAT, vv, local_n_bodies*3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(m, local_n_bodies, MPI_FLOAT, mm, local_n_bodies, MPI_FLOAT, 0, MPI_COMM_WORLD);

  if (my_rank == 0) {
    for (uint64_t i = 0; i < n_bodies*3; i+=3) {
      std::cout << vv[i] << std::endl;
      std::cout << vv[i + 1] << std::endl;
      std::cout << vv[i + 2] << std::endl;
    }

    for (uint64_t i = 0; i < n_bodies*3; i+=3) {
      std::cout << vv[i] << std::endl;
      std::cout << vv[i + 1] << std::endl;
      std::cout << vv[i + 2] << std::endl;
    }

    for (uint64_t i = 0; i < n_bodies; ++i) {
      std::cout << mm[i] << std::endl;
    }
  }

  MPI_Finalize();
}

int main(int argc, char **argv) {
  mpi_ap(12, 1.0f);
}