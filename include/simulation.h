#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <memory>
#include <random>

#include "data_format.h"
#include "nbody.h"
#include "sequential_ap.h"

// Implement the AllPairsSimulation subclass
class Simulation {
 public:
  Simulation() {
    randomizeBodies(this->simulation);
  }

 Simulation(const std::string& fname) {

   float *masses = nullptr, *positions = nullptr, *velocities = nullptr;
   float grav_const;
   auto n_bodies= LoadFromCSVConfiguration(fname, &masses, &positions, &velocities, grav_const);
   this->simulation = SequentialAP(n_bodies, masses, positions, velocities, grav_const);
  }

  void start(float time, float dt, const std::string& fname) {
    for (float t = 0.0; t < time; t += dt) {
      if (fname.length())
        ExportToCSV(this->simulation.positions, this->simulation.n_bodies, fname);
      this->simulation.Update(dt);
    }
  }

  void randomizeBodies(SequentialAP& nBody) {
      auto scale = std::max<float>(1.0f, nBody.n_bodies / (1024.0f));

      int i = 0;

      std::mt19937 mt;
      mt.seed(0);
      std::uniform_real_distribution<float> dt(-1,1);

      while (i < nBody.n_bodies * 3) {
        nBody.masses[i/3] = 1.0f;

        nBody.positions[i]   = dt(mt);
        nBody.positions[i+1] = dt(mt);
        nBody.positions[i+2] = dt(mt);

        auto lenSqr = nBody.positions[i] * nBody.positions[i]
            + nBody.positions[i+1] * nBody.positions[i+1]
            + nBody.positions[i+2] * nBody.positions[i+2];

        if (lenSqr > 1)
          continue;

        nBody.velocities[i]   = dt(mt);;
        nBody.velocities[i+1] = dt(mt);;
        nBody.velocities[i+2] = dt(mt);;

        lenSqr = nBody.velocities[i] * nBody.velocities[i]
            + nBody.velocities[i+1] * nBody.velocities[i+1]
            + nBody.velocities[i+2] * nBody.velocities[i+2];

        if (lenSqr > 1)
          continue;

        i+=3;
      }
  }


 private:
  SequentialAP simulation;
  //std::unique_ptr<NBody> simulation;
};

// Implement the AllPairsSimulation subclass
class ParallelAllPairsSimulation : public Simulation {
  // public:
  //  void Update() override {

  // }
};

#endif
