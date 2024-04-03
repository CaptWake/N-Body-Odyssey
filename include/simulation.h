#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <memory>
#include <random>

#include "data_format.h"
#include "nbody.h"
#include "sequential_ap.h"
#include "sequential_ap_avx.h"

class Simulation {
 public:
  Simulation() = default;

  Simulation(const std::string& mode, const std::string& fname) {
   if (mode == "SEQ") {
     this->sim_ap_seq = SequentialAP(fname);
   } else if (mode == "SEQ_AVX") {
     this->sim_ap_avx_seq = SequentialAPAVX(fname);
   }
   this->mode = mode;
  }

  void start(float time, float dt, const std::string& fname="") {
    for (float t = 0.0; t < time; t += dt) {
      if (this->mode == "SEQ") {
        if (!fname.empty())
          this->sim_ap_seq.LogsToCSV(fname);
        this->sim_ap_seq.Update(dt);
      }
      else if (this->mode == "SEQ_AVX") {
        if (!fname.empty())
          this->sim_ap_avx_seq.LogsToCSV(fname);
        this->sim_ap_avx_seq.Update(dt);
      }
    }
  }

  /*
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
  */
 protected:
  SequentialAP sim_ap_seq;
  SequentialAPAVX sim_ap_avx_seq;
  std::string mode;
};

// Implement the AllPairsSimulation subclass
class ParallelAllPairsSimulation : public Simulation {
  // public:
  //  void Update() override {

  // }
};

#endif
