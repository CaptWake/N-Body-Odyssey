#ifndef MPI_AP_H_
#define MPI_AP_H_

template <typename T>
void MPIAPSimulate(int n, T dt, T tEnd);
template <typename T>
void MPIAPUpdate(int localN, int n, const T* __restrict__ m,
                 const T* __restrict__ p, T* a, T* __restrict__ m_rec,
                 T* __restrict__ p_rec);

#endif
