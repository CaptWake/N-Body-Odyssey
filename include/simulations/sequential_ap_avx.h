#ifndef SEQUENTIAL_AP_AVX_H_
#define SEQUENTIAL_AP_AVX_H_

#include <cstdint>

// Update //
void SequentialAPAVXUpdate(float dt);

// Simulate laepfrog//
void SequentialAPAVXSimulate(uint64_t n, float dt, float tEnd, uint64_t seed);

#endif
