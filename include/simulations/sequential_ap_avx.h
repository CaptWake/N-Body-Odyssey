#ifndef SEQUENTIAL_AP_AVX_H_
#define SEQUENTIAL_AP_AVX_H_

#include <cstdint>

// Update //
#ifdef FLOAT
void SequentialAPAVXUpdate(float dt);
#else
void SequentialAPAVXUpdate(double dt);
#endif

// Simulate laepfrog//
template<typename T>
void SequentialAPAVXSimulate(uint64_t n, T dt, T tEnd, uint64_t seed);

#endif
