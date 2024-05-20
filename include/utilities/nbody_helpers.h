//
// Created by maste on 5/4/2024.
//

#ifndef NBODY_HELPERS_H_
#define NBODY_HELPERS_H_

#ifdef DOUBLE
#define MPI_TYPE MPI_DOUBLE
#define MY_T double
#else
#define MPI_TYPE MPI_FLOAT
#define MY_T float
#endif

#include <stdlib.h>

#include <cstdint>
#include <random>
#include <string>
#ifdef OMP
#include <omp.h>
#endif

constexpr float _G = 1;
constexpr float _M = 1;
constexpr float _SOFTENING = 0.025f;

/* pi */
static const double _PI = 2.0 * asin(1);

template <typename T>
static inline float fdrand() {
  return ((T)rand()) / (T)RAND_MAX;
}

// Function to export bodies information to CSV file

// void LogsToCSV(const std::string& filename);

// Initialization taken from https://github.com/alexgbrandt/Parallel-NBody/

template <typename T>
T inline Ep(uint64_t n, T *m, T *p) {
  T Epot = 0.0f;
  T D, x, y, z;
  uint64_t i, j;
  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
      x = p[3 * i + 0] - p[3 * j + 0];
      y = p[3 * i + 1] - p[3 * j + 1];
      z = p[3 * i + 2] - p[3 * j + 2];
      D = sqrtf(x * x + y * y + z * z);
      Epot += -1.0f * m[i] * m[j] / D;
    }
  }
  return Epot;
}

template <typename T>
T inline Ek(uint64_t n, T *m, T *v) {
  T Ekin = 0.0;
  uint64_t i;
  for (i = 0; i < n; ++i) {
    Ekin += 0.5f * m[i] *
            (v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] +
             v[3 * i + 2] * v[3 * i + 2]);
  }
  return Ekin;
}

template <typename T>
static inline void scale3NArray(uint64_t n, T *m, T scale) {
  for (uint64_t i = 0; i < 3 * n; ++i) {
    m[i] *= scale;
  }
}

template <typename T>
void InitMassU(uint64_t n, T *m) {
  T mi = _M / n;
  uint64_t i;

  for (i = 0; i < n; ++i) {
    m[i] = mi;
  }
}

template <typename T>
void InitPosU(uint64_t n, T *p) {
  T R, X, Y;
  uint64_t i;
  for (i = 0; i < n; ++i) {
    R = fdrand<T>();
    X = acos(1.0 - 2.0 * fdrand<T>());
    Y = fdrand<T>() * 2.0 * _PI;

    // https://www.researchgate.net/figure/Figure-A1-Spherical-coordinates_fig8_284609648
    p[3 * i + 0] = R * sin(X) * cos(Y);
    p[3 * i + 1] = R * sin(X) * sin(Y);
    p[3 * i + 2] = R * cos(X);
  }
}

template <typename T>
void InitVelU(uint64_t n, T *v) {
  uint64_t i;
  for (i = 0; i < n; ++i) {
    v[3 * i] = (1.0f - 2.0 * fdrand<T>());
    v[3 * i + 1] = (1.0f - 2.0 * fdrand<T>());
    v[3 * i + 2] = (1.0f - 2.0 * fdrand<T>());
  }
}

template <typename T>
void InitAccU(uint64_t n, T *a) {
  uint64_t i;
  for (i = 0; i < n; ++i) {
    a[3 * i] = 0.0;
    a[3 * i + 1] = 0.0;
    a[3 * i + 2] = 0.0;
  }
}

template <typename T>
void Move2Center(uint64_t n, T *m, T *p, T *v) {
  T px = 0.0, py = 0.0, pz = 0.0;
  T vx = 0.0, vy = 0.0, vz = 0.0;
  T mi;
  uint64_t i;

  for (i = 0; i < n; ++i) {
    mi = m[i];
    px += p[3 * i] * mi;
    py += p[3 * i + 1] * mi;
    pz += p[3 * i + 2] * mi;

    vx += v[3 * i] * mi;
    vy += v[3 * i + 1] * mi;
    vz += v[3 * i + 2] * mi;
  }

  px /= _M;
  py /= _M;
  pz /= _M;
  vx /= _M;
  vy /= _M;
  vz /= _M;

  for (i = 0; i < n; ++i) {
    p[3 * i] -= px;
    p[3 * i + 1] -= py;
    p[3 * i + 2] -= pz;
    v[3 * i] -= vx;
    v[3 * i + 1] -= vy;
    v[3 * i + 2] -= vz;
  }
}

template <typename T>
void RescaleEnergy(uint64_t n, T *m, T *p, T *v) {
  // Aarseth, 2003, Algorithm 7.2.
  T Epot = Ep<T>(n, m, p);
  T Ekin = Ek<T>(n, m, v);
  T virialRatio = 0.5;
  T Qv = sqrt(virialRatio * fabs(Epot) / Ekin);
  scale3NArray<T>(n, v, Qv);
  T beta = fabs((1 - virialRatio) * Epot / (Epot + Ekin));

  scale3NArray<T>(n, p, beta);
  scale3NArray<T>(n, v, 1.0f / (sqrt(beta)));

  // After first scale Ekin is -0.5Epot but E0 != -0.25.
  // So just scale up or down as needed.
  Epot = Ep<T>(n, m, p);
  beta = Epot / -0.5f;
  scale3NArray<T>(n, p, beta);
  scale3NArray<T>(n, v, 1.0f / sqrt(beta));
}

template <typename T>
void InitAos(const uint64_t n, T *m, T *p, T *v, T *a = nullptr) {
  // Initialize masses equally
  InitMassU<T>(n, m);

  // Initialize position with uniform distribution
  InitPosU(n, p);

  // Initialize velocities with uniform distribution
  InitVelU(n, v);

  if (a != nullptr)
    // Initialize masses equally
    InitAccU<T>(n, a);

  // Translate bodies to move the center of mass on center of the coordinate
  // system
  Move2Center<T>(n, m, p, v);

  // Rescale energy
  RescaleEnergy<T>(n, m, p, v);
}

// SOA

template <typename T>
static inline void scaleArray(uint64_t n, T *m, T scale) {
  for (uint64_t i = 0; i < n; ++i) {
    m[i] *= scale;
  }
}

template <typename T>
float inline EpSoa(uint64_t n, const T *m, const T *px, const T *py,
                   const T *pz) {
  T Epot = 0.0;
  T D, x, y, z;
  uint64_t i, j;
  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
      x = px[i] - px[j];
      y = py[i] - py[j];
      z = pz[i] - pz[j];
      D = sqrt(x * x + y * y + z * z);
      Epot += -1.0 * m[i] * m[j] / D;
    }
  }
  return Epot;
}

template <typename T>
float inline EkSoa(uint64_t n, const T *m, const T *vx, const T *vy,
                   const T *vz) {
  T Ekin = 0.0;
  uint64_t i;
  for (i = 0; i < n; ++i) {
    Ekin += 0.5 * m[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
  }
  return Ekin;
}

template <typename T>
void InitPosUSoa(uint64_t n, T *px, T *py, T *pz) {
  T R, X, Y;
  uint64_t i;
  for (i = 0; i < n; ++i) {
    R = fdrand<T>();
    X = acos(1.0 - 2.0 * fdrand<T>());
    Y = fdrand<T>() * 2.0 * _PI;

    // https://www.researchgate.net/figure/Figure-A1-Spherical-coordinates_fig8_284609648
    px[i] = R * sin(X) * cos(Y);
    py[i] = R * sin(X) * sin(Y);
    pz[i] = R * cos(X);
  }
}

template <typename T>
void InitVelUSoa(uint64_t n, T *vx, T *vy, T *vz) {
  uint64_t i;
  for (i = 0; i < n; ++i) {
    vx[i] = (1.0 - 2.0 * fdrand<T>());
    vy[i] = (1.0 - 2.0 * fdrand<T>());
    vz[i] = (1.0 - 2.0 * fdrand<T>());
  }
}

template <typename T>
void InitAccUSoa(uint64_t n, T *ax, T *ay, T *az) {
  uint64_t i;
  for (i = 0; i < n; ++i) {
    ax[i] = 0.0;
    ay[i] = 0.0;
    az[i] = 0.0;
  }
}

template <typename T>
void Move2CenterSoa(uint64_t n, T *m, T *px, T *py, T *pz, T *vx, T *vy,
                    T *vz) {
  T ppx = 0.0, ppy = 0.0, ppz = 0.0;
  T vvx = 0.0, vvy = 0.0, vvz = 0.0;
  T mi;
  uint64_t i;

  for (i = 0; i < n; ++i) {
    mi = m[i];
    ppx += px[i] * mi;
    ppy += py[i] * mi;
    ppz += pz[i] * mi;

    vvx += vx[i] * mi;
    vvy += vy[i] * mi;
    vvz += vz[i] * mi;
  }

  ppx /= _M;
  ppy /= _M;
  ppz /= _M;
  vvx /= _M;
  vvy /= _M;
  vvz /= _M;

  for (i = 0; i < n; ++i) {
    px[i] -= ppx;
    py[i] -= ppy;
    pz[i] -= ppz;
    vx[i] -= vvx;
    vy[i] -= vvy;
    vz[i] -= vvz;
  }
}

template <typename T>
void RescaleEnergySoa(uint64_t n, T *m, T *px, T *py, T *pz, T *vx, T *vy,
                      T *vz) {
  // Aarseth, 2003, Algorithm 7.2.
  T Epot = EpSoa<T>(n, m, px, py, pz);
  T Ekin = EkSoa<T>(n, m, vx, vy, vz);
  T virialRatio = 0.5;
  T Qv = sqrt(virialRatio * fabs(Epot) / Ekin);
  scaleArray<T>(n, vx, Qv);
  scaleArray<T>(n, vy, Qv);
  scaleArray<T>(n, vz, Qv);
  T beta = fabs((1 - virialRatio) * Epot / (Epot + Ekin));

  scaleArray<T>(n, px, beta);
  scaleArray<T>(n, py, beta);
  scaleArray<T>(n, pz, beta);
  scaleArray<T>(n, vx, 1.0 / (sqrt(beta)));
  scaleArray<T>(n, vy, 1.0 / (sqrt(beta)));
  scaleArray<T>(n, vz, 1.0 / (sqrt(beta)));

  // After first scale Ekin is -0.5Epot but E0 != -0.25.
  // So just scale up or down as needed.
  Epot = EpSoa<T>(n, m, px, py, pz);
  beta = Epot / -0.5;
  scaleArray<T>(n, px, beta);
  scaleArray<T>(n, py, beta);
  scaleArray<T>(n, pz, beta);
  scaleArray<T>(n, vx, 1.0 / sqrt(beta));
  scaleArray<T>(n, vy, 1.0 / sqrt(beta));
  scaleArray<T>(n, vz, 1.0 / sqrt(beta));
}

template <typename T>
void InitSoa(const uint64_t n, T *m, T *px, T *py, T *pz, T *vx, T *vy, T *vz) {
  // Initialize masses equally
  InitMassU<T>(n, m);

  // Initialize position with uniform distribution
  InitPosUSoa<T>(n, px, py, pz);

  // Initialize velocities with uniform distribution
  InitVelUSoa<T>(n, vx, vy, vz);

  // Translate bodies to move the center of mass on center of the coordinate
  // system
  Move2CenterSoa<T>(n, m, px, py, pz, vx, vy, vz);

  // Rescale energy
  RescaleEnergySoa<T>(n, m, px, py, pz, vx, vy, vz);
}

#ifdef OMP
// OMP
void SetScheduleType(const std::string &schedule_type, int chunk_size) {
  if (schedule_type == "static") {
    omp_set_schedule(omp_sched_static, chunk_size);
  } else {
    omp_set_schedule(omp_sched_dynamic, chunk_size);
  }
}

void SetNumThread(int num_threads) { omp_set_num_threads(num_threads); }
#endif
#endif  // NBODY_HELPERS_H_
