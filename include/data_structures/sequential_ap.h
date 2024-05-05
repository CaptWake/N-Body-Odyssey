#ifndef SEQUENTIALAP_H_
#define SEQUENTIALAP_H_

#include <cstdint>
#include <iostream>
#include <random>


  // Initialization //
  void InitSequentialAP(uint64_t n, float *m, float *p, float *v);

  // Update //
  void SequentialAPUpdate(float dt);

  // Simulate laepfrog//
  void SequentialAPSimulate(uint64_t n, float dt, float tEnd, uint64_t seed);

  // Simulate kick drift kick //
  void SequentialAPSimulate2(uint64_t n, float dt, float tEnd, uint64_t seed);



#endif
