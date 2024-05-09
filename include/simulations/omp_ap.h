#ifndef OMP_AP_H_
#define OMP_AP_H_

#include <cstdint>
#include <random>
#include <string>
#include <vector>

// Update //
void OMPAPUpdate(float dt);

// Simulate laepfrog//
void OMPAPSimulate(uint64_t n, float dt, float tEnd, uint64_t seed);

// Simulate kick drift kick //
void OMPAPSimulate2(uint64_t n, float dt, float tEnd, uint64_t seed);

#endif
