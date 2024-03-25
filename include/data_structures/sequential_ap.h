//
// Created by maste on 3/21/2024.
//

#ifndef SEQUENTIALAP_H_
#define SEQUENTIALAP_H_

#include <cstdint>
#include <iostream>
#include <cmath>

#include "nbody.h"

class SequentialAP : public NBody {
  public:
   SequentialAP() { G=1; masses=nullptr; velocities=nullptr; positions=nullptr; n_bodies=0; }
   SequentialAP(uint64_t n_bodies, float grav_const) {
     this->n_bodies = n_bodies;
     this->masses = new float[n_bodies];
     this->positions = new float[n_bodies * 3];
     this->velocities = new float[n_bodies * 3];
     this->G = grav_const;
   }

   // Update//
   void Update(float dt) override;

   friend std::ostream& operator<<(std::ostream& os, const SequentialAP& nbody);

  ~SequentialAP() {
    delete[] masses;
    delete[] positions;
    delete[] velocities;
  }

  float *masses;
  float *positions;
  float *velocities;
  uint64_t n_bodies;
  // mass, position, velocity, <pad>, mass, position, velocity, <pad>
  float G;
};

#endif
