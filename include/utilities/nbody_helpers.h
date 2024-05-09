//
// Created by maste on 5/4/2024.
//

#ifndef NBODY_HELPERS_H_
#define NBODY_HELPERS_H_

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

static inline float frand() { return ((float)rand()) / (float)RAND_MAX; }

// Function to export bodies information to CSV file

// void LogsToCSV(const std::string& filename);

// Initialization taken from https://github.com/alexgbrandt/Parallel-NBody/

float inline Ep(uint64_t n, float *m, float *p) {
  float Epot = 0.0f;
  float D, x, y, z;
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

float inline Ek(uint64_t n, float *m, float *v) {
  float Ekin = 0.0;
  uint64_t i;
  for (i = 0; i < n; ++i) {
    Ekin += 0.5f * m[i] *
            (v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] +
             v[3 * i + 2] * v[3 * i + 2]);
  }
  return Ekin;
}

static inline void scale3NArray(uint64_t n, float *m, float scale) {
  for (uint64_t i = 0; i < 3 * n; ++i) {
    m[i] *= scale;
  }
}

void InitMassU(uint64_t n, float *m) {
  float mi = _M / n;
  uint64_t i;

  for (i = 0; i < n; ++i) {
    m[i] = mi;
  }
}

void InitPosU(uint64_t n, float *p) {
  float R, X, Y;
  uint64_t i;
  for (i = 0; i < n; ++i) {
    R = frand();
    X = acosf(1.0f - 2.0f * frand());
    Y = frand() * 2.0f * _PI;

    // https://www.researchgate.net/figure/Figure-A1-Spherical-coordinates_fig8_284609648
    p[3 * i + 0] = R * sinf(X) * cosf(Y);
    p[3 * i + 1] = R * sinf(X) * sinf(Y);
    p[3 * i + 2] = R * cosf(X);
  }
}

void InitVelU(uint64_t n, float *v) {
  uint64_t i;
  for (i = 0; i < n; ++i) {
    v[3 * i] = (1.0f - 2.0f * frand());
    v[3 * i + 1] = (1.0f - 2.0f * frand());
    v[3 * i + 2] = (1.0f - 2.0f * frand());
  }
}

void InitAccU(uint64_t n, float *a) {
  uint64_t i;
  for (i = 0; i < n; ++i) {
    a[3 * i] = 0.0f;
    a[3 * i + 1] = 0.0f;
    a[3 * i + 2] = 0.0f;
  }
}

void Move2Center(uint64_t n, float *m, float *p, float *v) {
  float px = 0.0f, py = 0.0f, pz = 0.0f;
  float vx = 0.0f, vy = 0.0f, vz = 0.0f;
  float mi;
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

void RescaleEnergy(uint64_t n, float *m, float *p, float *v) {
  // Aarseth, 2003, Algorithm 7.2.
  float Epot = Ep(n, m, p);
  float Ekin = Ek(n, m, v);
  float virialRatio = 0.5f;
  float Qv = sqrtf(virialRatio * fabsf(Epot) / Ekin);
  scale3NArray(n, v, Qv);
  float beta = fabsf((1 - virialRatio) * Epot / (Epot + Ekin));

  scale3NArray(n, p, beta);
  scale3NArray(n, v, 1.0f / (sqrtf(beta)));

  // After first scale Ekin is -0.5Epot but E0 != -0.25.
  // So just scale up or down as needed.
  Epot = Ep(n, m, p);
  beta = Epot / -0.5f;
  scale3NArray(n, p, beta);
  scale3NArray(n, v, 1.0f / sqrtf(beta));
}

void InitAos(const uint64_t n, float *m, float *p, float *v,
             float *a = nullptr) {
  // Initialize masses equally
  InitMassU(n, m);

  // Initialize position with uniform distribution
  InitPosU(n, p);

  // Initialize velocities with uniform distribution
  InitVelU(n, v);

  if (a != nullptr)
    // Initialize masses equally
    InitAccU(n, a);

  // Translate bodies to move the center of mass on center of the coordinate
  // system
  Move2Center(n, m, p, v);

  // Rescale energy
  RescaleEnergy(n, m, p, v);
}

// SOA

static inline void scaleArray(uint64_t n, float *m, float scale) {
  for (uint64_t i = 0; i < n; ++i) {
    m[i] *= scale;
  }
}

float inline EpSoa(uint64_t n, const float *m, const float *px, const float *py,
                   const float *pz) {
  float Epot = 0.0f;
  float D, x, y, z;
  uint64_t i, j;
  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
      x = px[i] - px[j];
      y = py[i] - py[j];
      z = pz[i] - pz[j];
      D = sqrtf(x * x + y * y + z * z);
      Epot += -1.0f * m[i] * m[j] / D;
    }
  }
  return Epot;
}

float inline EkSoa(uint64_t n, const float *m, const float *vx, const float *vy,
                   const float *vz) {
  float Ekin = 0.0;
  uint64_t i;
  for (i = 0; i < n; ++i) {
    Ekin += 0.5f * m[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
  }
  return Ekin;
}

void InitPosUSoa(uint64_t n, float *px, float *py, float *pz) {
  float R, X, Y;
  uint64_t i;
  for (i = 0; i < n; ++i) {
    R = frand();
    X = acosf(1.0f - 2.0f * frand());
    Y = frand() * 2.0f * _PI;

    // https://www.researchgate.net/figure/Figure-A1-Spherical-coordinates_fig8_284609648
    px[i] = R * sinf(X) * cosf(Y);
    py[i] = R * sinf(X) * sinf(Y);
    pz[i] = R * cosf(X);
  }
}

void InitVelUSoa(uint64_t n, float *vx, float *vy, float *vz) {
  uint64_t i;
  for (i = 0; i < n; ++i) {
    vx[i] = (1.0f - 2.0f * frand());
    vy[i] = (1.0f - 2.0f * frand());
    vz[i] = (1.0f - 2.0f * frand());
  }
}

void InitAccUSoa(uint64_t n, float *ax, float *ay, float *az) {
  uint64_t i;
  for (i = 0; i < n; ++i) {
    ax[i] = 0.0f;
    ay[i] = 0.0f;
    az[i] = 0.0f;
  }
}

void Move2CenterSoa(uint64_t n, float *m, float *px, float *py, float *pz,
                    float *vx, float *vy, float *vz) {
  float ppx = 0.0f, ppy = 0.0f, ppz = 0.0f;
  float vvx = 0.0f, vvy = 0.0f, vvz = 0.0f;
  float mi;
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

void RescaleEnergySoa(uint64_t n, float *m, float *px, float *py, float *pz,
                      float *vx, float *vy, float *vz) {
  // Aarseth, 2003, Algorithm 7.2.
  float Epot = EpSoa(n, m, px, py, pz);
  float Ekin = EkSoa(n, m, vx, vy, vz);
  float virialRatio = 0.5f;
  float Qv = sqrtf(virialRatio * fabsf(Epot) / Ekin);
  scaleArray(n, vx, Qv);
  scaleArray(n, vy, Qv);
  scaleArray(n, vz, Qv);
  float beta = fabsf((1 - virialRatio) * Epot / (Epot + Ekin));

  scaleArray(n, px, beta);
  scaleArray(n, py, beta);
  scaleArray(n, pz, beta);
  scaleArray(n, vx, 1.0f / (sqrtf(beta)));
  scaleArray(n, vy, 1.0f / (sqrtf(beta)));
  scaleArray(n, vz, 1.0f / (sqrtf(beta)));

  // After first scale Ekin is -0.5Epot but E0 != -0.25.
  // So just scale up or down as needed.
  Epot = EpSoa(n, m, px, py, pz);
  beta = Epot / -0.5f;
  scaleArray(n, px, beta);
  scaleArray(n, py, beta);
  scaleArray(n, pz, beta);
  scaleArray(n, vx, 1.0f / sqrtf(beta));
  scaleArray(n, vy, 1.0f / sqrtf(beta));
  scaleArray(n, vz, 1.0f / sqrtf(beta));
}

void InitSoa(const uint64_t n, float *m, float *px, float *py, float *pz,
             float *vx, float *vy, float *vz) {
  // Initialize masses equally
  InitMassU(n, m);

  // Initialize position with uniform distribution
  InitPosUSoa(n, px, py, pz);

  // Initialize velocities with uniform distribution
  InitVelUSoa(n, vx, vy, vz);

  // Translate bodies to move the center of mass on center of the coordinate
  // system
  Move2CenterSoa(n, m, px, py, pz, vx, vy, vz);

  // Rescale energy
  RescaleEnergySoa(n, m, px, py, pz, vx, vy, vz);
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
