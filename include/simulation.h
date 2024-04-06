#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <memory>
#include <random>

#include "data_format.h"
#include "nbody.h"
#include "sequential_ap.h"
#include "sequential_ap_avx.h"
#include "sequential_bh.h"

class Simulation {
 public:
  Simulation() = default;

  Simulation(const std::string& algorithm, const std::string& mode, const std::string& fname) {
   if (algorithm == "AP") {
      if (mode == "SEQ") {
       this->sim_ap_seq = SequentialAP(2048, 1.0f); // SequentialAP(fname);
     } else if (mode == "SEQ_AVX") {
       this->sim_ap_avx_seq = SequentialAPAVX(2048, 1.0f); // SequentialAPAVX(fname);
     }
   // Barnes Hut
   } else {
     if (mode == "SEQ") {
       this->sim_bh_seq = SequentialBH(2048, 1.0f, 0.9f);
     } else if (mode == "SEQ_AVX") {
       this->sim_ap_avx_seq = SequentialAPAVX(2048, 1.0f); // SequentialAPAVX(fname);
     }
   }
   this->algorithm = algorithm;
   this->mode = mode;
  }


  void start(float time, float dt, const std::string& fname="") {
    float t = 0.0;
    while (t < time) {
      if (this->algorithm == "AP") {
        if (this->mode == "SEQ") {
          if (!fname.empty())
            this->sim_ap_seq.LogsToCSV(fname);
          this->sim_ap_seq.Update(dt);
        } else if (this->mode == "SEQ_AVX") {
          if (!fname.empty())
            this->sim_ap_avx_seq.LogsToCSV(fname);
          this->sim_ap_avx_seq.Update(dt);
        }
      // assuming Barnes Hut
      } else {
        if (this->mode == "SEQ") {
          if (!fname.empty())
            this->sim_bh_seq.LogsToCSV(fname);
          this->sim_bh_seq.Update(dt);
        } else if (this->mode == "SEQ_AVX") {
          if (!fname.empty())
            this->sim_ap_avx_seq.LogsToCSV(fname);
          this->sim_ap_avx_seq.Update(dt);
        }
      }
      t += dt;
    }
  }

 protected:
  SequentialAP sim_ap_seq;
  SequentialAPAVX sim_ap_avx_seq;
  SequentialBH sim_bh_seq;
  std::string algorithm, mode;
};

#endif
