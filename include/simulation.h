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

    //particles.emplace_back(1, vec3(-0.97000436, 0.24308753, 0), vec3(0.4662036850, 0.4323657300,0));
    //particles.emplace_back(1, vec3(-0, 0, 0), vec3(-0.93240737, -0.86473146,0));
    //particles.emplace_back(1, vec3(0.97000436, -0.24308753, 0), vec3(0.4662036850, 0.4323657300,0));

    float *masses = nullptr, *positions = nullptr, *velocities = nullptr;
    float grav_const;
    auto n_bodies= LoadFromCSVConfiguration("/home/ste/Downloads/SolarSystem.csv", &masses, &positions, &velocities, grav_const);
    this->simulation = SequentialAP(n_bodies, masses, positions, velocities, grav_const);
/*
    this->simulation = SequentialAP(3, 1.0f);
    this->simulation.masses[0] = 1;

    this->simulation.positions[0] = -0.97000436;
    this->simulation.positions[1] = 0.24308753;
    this->simulation.positions[2] = 0;

    this->simulation.velocities[0] = 0.4662036850;
    this->simulation.velocities[1] = 0.4323657300;
    this->simulation.velocities[2] = 0;

    this->simulation.masses[1] = 1;

    this->simulation.positions[3] = -0;
    this->simulation.positions[4] = 0;
    this->simulation.positions[5] = 0;

    this->simulation.velocities[3] = -0.93240737;
    this->simulation.velocities[4] = -0.86473146;
    this->simulation.velocities[5] = 0;

    this->simulation.masses[2] = 1;

    this->simulation.positions[6] = 0.97000436;
    this->simulation.positions[7] = -0.24308753;
    this->simulation.positions[8] = 0;

    this->simulation.velocities[6] = 0.4662036850;
    this->simulation.velocities[7] = 0.4323657300;
    this->simulation.velocities[8] = 0;

    // randomizeBodies(this->simulation);
*/
  }

  void start(float time, float dt) {
    for (float t = 0.0; t < time; t += dt) {
      ExportToCSV(this->simulation.positions, this->simulation.n_bodies, "/home/ste/Documents/SCPD-Project/src/results_correct.csv");
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
