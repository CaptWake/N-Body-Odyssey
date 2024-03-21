#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "data_format.h"
#include "nbody.h"

// Define the Simulation interface
class Simulation {
 public:
  virtual void start(double time, double dt) = 0;  // Pure virtual method
  virtual ~Simulation() {}  // Virtual destructor for polymorphic behavior
};

// Implement the AllPairsSimulation subclass
class SequentialAllPairsSimulation : public Simulation {
 public:
  SequentialAllPairsSimulation() {
    std::vector<Particle> particles;
    particles.emplace_back(1, vec3(-0.97000436, 0.24308753, 0), vec3(0.4662036850, 0.4323657300,0));
    particles.emplace_back(1, vec3(-0, 0, 0), vec3(-0.93240737, -0.86473146,0));
    particles.emplace_back(1, vec3(0.97000436, -0.24308753, 0), vec3(0.4662036850, 0.4323657300,0));
    this->simulation = NBody(particles, 1);
  }
  void start(double time, double dt) override {
    int i = 0;
    for (double t = 0.0; t < time; t += dt) {
      std::cout << "Iteration " << i << std::endl;
      ExportToCSV(this->simulation.GetParticles(),
                  "/home/ste/Documents/SCPD-Project/src/results.csv");
      this->simulation.update(t);
      i += 1;
    }
  }

 private:
  NBody simulation;
};

// Implement the AllPairsSimulation subclass
class ParallelAllPairsSimulation : public Simulation {
  // public:
  //  void update() override {

  // }
};

#endif
