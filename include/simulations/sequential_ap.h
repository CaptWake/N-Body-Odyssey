#ifndef SEQUENTIALAP_H_
#define SEQUENTIALAP_H_

#include <cstdint>
#include <iostream>
#include <random>

// Update //
template <typename T>
void SequentialAPUpdate(T dt);

// Simulate laepfrog//
template <typename T>
void SequentialAPSimulate(uint64_t n, T dt, T tEnd, uint64_t seed);

// Simulate kick drift kick //
template <typename T>
void SequentialAPSimulate2(uint64_t n, T dt, T tEnd, uint64_t seed);

#endif
