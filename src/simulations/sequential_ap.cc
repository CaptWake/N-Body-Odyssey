#include "sequential_ap.h"

#include <cmath>
#include <fstream>
#include "nbody_helpers.h"
#include "time_utils.h"

void InitSequentialAP(const uint64_t n, float *m, float *p, float *v, float *a){

  // Initialize masses equally
  InitMassU(n, m);

  // Initialize position with uniform distribution
  InitPosU(n, p);

  // Initialize velocities with uniform distribution
  InitVelU(n, v);

  // Initialize masses equally
  InitAccU(n, a);

  // Translate bodies to move the center of mass on center of the coordinate system
  Move2Center(n, m, p, v);

  // Rescale energy
  RescaleEnergy(n, m, p, v);
}

void computeForce(uint64_t n, uint64_t i,
                  const float* m,
                  const float* p,
                  float* a){
  //         a[k] = -p[k] / (r2*sqrtr2);
  float* ai = a + 3*i;
  ai[0] = 0.0f; ai[1] = 0.0f; ai[2] = 0.0f;
  float px = p[3*i + 0];
  float py = p[3*i + 1];
  float pz = p[3*i + 2];
  float dx, dy, dz, D;
  uint64_t j;
  for (j = 0; j < i; ++j) {
    //really dx is other way around, be this way we can avoid -1.0* later.
    dx = p[3*j] - px;
    dy = p[3*j + 1] - py;
    dz = p[3*j + 2] - pz;
    D = dx*dx + dy*dy + dz*dz;
    D += 1e-09*1e-09;
    D = 1.0f / (D*sqrtf(D));
    ai[0] += m[j]*dx*D;
    ai[1] += m[j]*dy*D;
    ai[2] += m[j]*dz*D;
  }
  for (j = i+1; j < n; ++j) {
    dx = p[3*j] - px;
    dy = p[3*j + 1] - py;
    dz = p[3*j + 2] - pz;
    D = dx*dx + dy*dy + dz*dz;
    D += 1e-09* 1e-09;
    D = 1.0f / (D*sqrtf(D));
    ai[0] += m[j]*dx*D;
    ai[1] += m[j]*dy*D;
    ai[2] += m[j]*dz*D;
  }
}

static inline void performNBodyHalfStepA(uint64_t n, float dt,
                                         float* p,
                                         float* v,
                                         const float* a,
                                         const float* m){
  for (uint64_t i = 0; i < n; ++i) {
    //kick, drift
    v[3*i + 0] += 0.5f * a[3*i + 0] * dt;
    p[3*i + 0] += v[3*i + 0] * dt;
    v[3*i + 1] += 0.5f * a[3*i + 1] * dt;
    p[3*i + 1] += v[3*i + 1] * dt;
    v[3*i + 2] += 0.5f * a[3*i + 2] * dt;
    p[3*i + 2] += v[3*i + 2] * dt;
  }
}

static inline void performNBodyHalfStepB(uint64_t n, float dt,
                                        const float* p,
                                        float* v,
                                        const float* a,
                                        const float* m)
{
  for (uint64_t i = 0; i < n; ++i) {
    //kick
    v[3*i + 0] += 0.5f * a[3*i + 0] * dt;
    v[3*i + 1] += 0.5f * a[3*i + 1] * dt;
    v[3*i + 2] += 0.5f * a[3*i + 2] * dt;
  }
}


// copyright NVIDIA
void SequentialAPUpdate(const uint64_t n, float *m, float *p , float *v, const float dt) {
  for (uint64_t i = 0; i < n * 3; i += 3) {
    float fx = 0.0f;
    float fy = 0.0f;
    float fz = 0.0f;
    for (uint64_t j = i + 1; j < n * 3; j += 3) {
      auto m2_id = j / 3;
      // compute distance pair
      auto dx = p[j] - p[i];

      auto dy = p[j + 1] - p[i + 1];
      auto dz = p[j + 2] - p[i + 2];

      auto d = dx * dx + dy * dy + dz * dz + 1e-9f * 1e-9f;
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

  for (uint64_t i = 0; i < n * 3; i += 3) {
    p[i] += v[i] * dt;
    p[i + 1] += v[i + 1] * dt;
    p[i + 2] += v[i + 2] * dt;
  }
}


void SequentialAPSimulate(uint64_t n, float dt, float tEnd, uint64_t seed){
  float *m = new float[n];
  float *p = new float[3 * n];
  float *v = new float[3 * n];
  float *a = new float[3 * n];

  // Init Bodies
  InitSequentialAP(n, m, p, v, a);

  // Simulation Loop
  for (float t = 0.0f; t < tEnd; t += dt){
    // Update Bodies
    SequentialAPUpdate(n, m, p, v, dt);
    float ek = Ek(n, m, v);
    float ep = Ep(n, m, p);
    std::cout << "Etot: " <<ek+ep <<std::endl;

    t += dt;
  }
}

void SequentialAPSimulate2(uint64_t n, float dt, float tEnd, uint64_t seed){
  float *m = new float[n];
  float *p = new float[3 * n];
  float *v = new float[3 * n];
  float *a = new float[3 * n];

  // Init Bodies
  InitSequentialAP(n, m, p, v, a);

  // Simulation Loop
  for (float t = 0.0f; t < tEnd; t += dt){
    // Update Bodies
    performNBodyHalfStepA(n, dt, p, v, a, m);
    for (uint64_t i = 0; i < n; ++i) {
      computeForce(n, i, m, p ,a);
    }
    performNBodyHalfStepB(n, dt, p, v, a, m);
    float ek = Ek(n, m, v);
    float ep = Ep(n, m, p);
    std::cout << "Etot: " <<ek+ep <<std::endl;

    t += dt;
  }
}

int main (int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  srand(1714995103);
  TIMERSTART(simulation)
  SequentialAPSimulate(std::stoul(argv[1]), 0.01, 10, 0);
  TIMERSTOP(simulation)
}