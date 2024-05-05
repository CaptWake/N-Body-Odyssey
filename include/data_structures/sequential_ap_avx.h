#ifndef SEQUENTIAL_AP_AVX_H_
#define SEQUENTIAL_AP_AVX_H_

#include <cstdint>
#include <iostream>
#include <random>

#include "fileIO.h"
// Initialization //
void InitSequentialAPAVX(uint64_t n, float *m, float *px, float *py, float*pz, float *vx, float *vy, float *vz);

// Update //
void SequentialAPAVXUpdate(float dt);

// Simulate laepfrog//
void SequentialAPAVXSimulate(uint64_t n, float dt, float tEnd, uint64_t seed);

#endif
