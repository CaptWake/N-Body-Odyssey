#ifndef SEQUENTIALAP_H_
#define SEQUENTIALAP_H_

#include <cstdint>
#include <iostream>
#include <random>
#include "nbody.h"

class SequentialAP : public NBody {
  public:
    SequentialAP() { G=1; masses=nullptr; velocities=nullptr; positions=nullptr; n_bodies=0; }

    explicit SequentialAP(const std::string& fname) {
      LoadFromCSVConfiguration(fname);
    }

    // generate random samples
    SequentialAP(uint64_t n_bodies, float grav_const) {
      static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
      // initialize Mersennes' twister using rd to generate the seed
      static std::mt19937 engine{0};//rd()};
      std::uniform_real_distribution<float> density(-1, 1);

      this->n_bodies   = n_bodies;
      this->masses     = new float[n_bodies];
      this->positions  = new float[n_bodies*3];
      this->velocities = new float[n_bodies*3];
      this->G = grav_const;

      for (uint64_t i = 0; i < n_bodies*3; i+=3) {
        this->masses[i/3]     = density(engine);

        this->positions[i]    = density(engine);
        this->positions[i+1]  = density(engine);
        this->positions[i+2]  = density(engine);

        this->velocities[i]   = density(engine);
        this->velocities[i+1] = density(engine);
        this->velocities[i+2] = density(engine);
      }
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
  void LogsToCSV(const std::string& filename) const;

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
