//
// Created by maste on 3/21/2024.
//

#ifndef SEQUENTIALAP_H_
#define SEQUENTIALAP_H_

#include <cstdint>
#include <iostream>
#include "nbody.h"

class SequentialAP : public NBody {
  public:
    SequentialAP() { G=1; masses=nullptr; velocities=nullptr; positions=nullptr; n_bodies=0; }

    explicit SequentialAP(const std::string& fname) {
      LoadFromCSVConfiguration(fname);
    }

   SequentialAP(uint64_t n_bodies, float grav_const) {
     this->n_bodies = n_bodies;
     this->masses = new float[n_bodies];
     this->positions = new float[n_bodies * 3];
     this->velocities = new float[n_bodies * 3];
     this->G = grav_const;
   }

  SequentialAP(uint64_t n_bodies, float* masses, float* positions, float* velocities, float grav_const) {
    this->n_bodies = n_bodies;
    this->masses = masses;
    this->positions = positions;
    this->velocities = velocities;
    this->G = grav_const;
  }

   //Move Constructor
   SequentialAP& operator=(SequentialAP&& old) noexcept {
       n_bodies=old.n_bodies;
       masses=old.masses;
       velocities=old.velocities;
       positions=old.positions;
       G=old.G;
       old.masses = nullptr;
       old.positions = nullptr;
       old.velocities = nullptr;
       return *this;
   }

  void LoadFromCSVConfiguration(const std::string& filename);

  // Function to export bodies information to CSV file
  void LogsToCSV(const std::string& filename);

   // Update//
   void Update(float dt) override;

   friend std::ostream& operator<<(std::ostream& os, const SequentialAP& nbody);

  ~SequentialAP() {
    free(masses);
    free(positions);
    free(velocities);
  }

  float *masses{};
  float *positions{};
  float *velocities{};
  uint64_t n_bodies{};
  // mass, position, velocity, <pad>, mass, position, velocity, <pad>
  float G{};
};

#endif
