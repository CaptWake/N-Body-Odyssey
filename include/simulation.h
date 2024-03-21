#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "data_format.h"
#include "nbody.h"


// Implement the AllPairsSimulation subclass
class Simulation {
 public:
  Simulation() {
    std::vector<Body> particles;
    particles.emplace_back(1, vec3(-0.97000436, 0.24308753, 0), vec3(0.4662036850, 0.4323657300,0));
    particles.emplace_back(1, vec3(-0, 0, 0), vec3(-0.93240737, -0.86473146,0));
    particles.emplace_back(1, vec3(0.97000436, -0.24308753, 0), vec3(0.4662036850, 0.4323657300,0));
    this->simulation = SequentialAPNBody(particles, 1);
  }
  void start(double time, double dt) {
    int i = 0;
    for (double t = 0.0; t < time; t += dt) {
      std::cout << "Iteration " << i << std::endl;
      ExportToCSV(this->simulation.GetBodies(),
                  "/home/ste/Documents/SCPD-Project/src/results.csv");
      this->simulation.update(t);
      i += 1;
    }
  }

 private:
  SequentialAPNBody simulation;
};

// Implement the AllPairsSimulation subclass
class ParallelAllPairsSimulation : public Simulation {
  // public:
  //  void update() override {

  // }
};

#endif
