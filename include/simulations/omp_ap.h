#ifndef OMP_AP_H_
#define OMP_AP_H_

#include <cstdint>
#include <random>
#include <string>
#include <vector>

template <typename T>
// Update //
void OMPAPUpdate(T dt);

template <typename T>
// Simulate laepfrog//
void OMPAPSimulate(int n, T dt, T tEnd);

template <typename T>
// Simulate kick drift kick //
void OMPAPSimulate2(int n, T dt, T tEnd);

#endif
