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
void SequentialAPSimulate(int n, T dt, T tEnd);

// Simulate kick drift kick //
template <typename T>
void SequentialAPSimulate2(int n, T dt, T tEnd);

#endif
