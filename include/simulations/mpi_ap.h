#ifndef MPI_AP_H_
#define MPI_AP_H_

  void MPIAPSimulate(int n, float dt, float tEnd, uint64_t seed);
  void MPIAPUpdate(int localN,
      int n,
      const float* __restrict__ m,
      const float* __restrict__ p,
      float* a,
      float* __restrict__ m_rec,
      float* __restrict__ p_rec);

#endif
