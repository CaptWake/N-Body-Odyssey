#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <memory>
#include <random>

#include "nbody.h"
#include "sequential_ap.h"
#include "sequential_ap_avx.h"
#include "sequential_bh.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

class Simulation {
 public:
  Simulation() = default;

  explicit Simulation(const json& sim_args) {
   if (sim_args["algorithm"] == "AP") {
      if (sim_args["mode"] == "SEQ") {
        if (sim_args.contains("ipath"))
          this->sim_ap_seq = SequentialAP(sim_args["ipath"]);
        else
          this->sim_ap_seq = SequentialAP(sim_args["n_bodies"], sim_args["grav_const"]);
     } else if (sim_args["mode"] == "SEQ_AVX") {
        if (sim_args.contains("ipath"))
          this->sim_ap_avx_seq = SequentialAPAVX(sim_args["ipath"]);
        else
          this->sim_ap_avx_seq = SequentialAPAVX(sim_args["n_bodies"], sim_args["grav_const"]);
     }
   // Barnes Hut
   } else {
     if (sim_args["mode"] == "SEQ") {
       if (sim_args.contains("ipath"))
         this->sim_bh_seq = SequentialBH(sim_args["ipath"]);
       else
         this->sim_bh_seq = SequentialBH(sim_args["n_bodies"], sim_args["grav_const"], sim_args["theta"]);
     }
   }
   this->algorithm = sim_args["algorithm"];
   this->mode = sim_args["mode"];
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
