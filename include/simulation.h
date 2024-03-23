#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <memory>

#include "data_format.h"
#include "nbody.h"
#include "sequential_ap.h"

// Implement the AllPairsSimulation subclass
class Simulation {
 public:
  Simulation() {
    std::vector<Body> particles;

    particles.emplace_back(1, vec3(-0.97000436, 0.24308753, 0), vec3(0.4662036850, 0.4323657300,0));
    particles.emplace_back(1, vec3(-0, 0, 0), vec3(-0.93240737, -0.86473146,0));
    particles.emplace_back(1, vec3(0.97000436, -0.24308753, 0), vec3(0.4662036850, 0.4323657300,0));

    //*3600/149597870700
    this->simulation = std::make_unique<SequentialAP>(SequentialAP(particles, 5.137597588561357e-07));
  }
  void start(double time, double dt) {
    for (double t = 0.0; t < time; t += dt) {
      this->simulation->Update(dt);
    }
  }

 private:
  std::unique_ptr<NBody> simulation;
};

// Implement the AllPairsSimulation subclass
class ParallelAllPairsSimulation : public Simulation {
  // public:
  //  void Update() override {

  // }
};

#endif
