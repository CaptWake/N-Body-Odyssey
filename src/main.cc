#include <argparse/argparse.hpp>
#include <fstream>

#include "octree.h"
#include "sequential_bh_avx.h"
#include "simulation.h"
#include "time_utils.h"

int main(int argc, char** argv) {
  //SequentialBHAVX(const uint64_t n_bodies, const float grav_const, float theta)
  auto sim_test = SequentialBHAVX(1024, 1, 0);
  auto sim_true = SequentialBH(1024, 1, 0);

  TIMERSTART(sim)
  for (float i = 0; i < 10; i+=0.01) {
    sim_true.Update(0.01);
  }
  TIMERSTOP(sim)

  TIMERSTART(avx)
  for (float i = 0; i < 10; i+=0.01) {
    sim_test.Update(0.01);
  }
  TIMERSTOP(avx)

  return 1;
  argparse::ArgumentParser program("N-Body-Simulator");
  program.add_argument("-c", "--config")
      .help("specify the configuration file path");

  try {
    program.parse_args(argc, argv);
  } catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  if (program.is_used("-c")) {
    auto config_path = program.get<std::string>("-c");
    std::ifstream f(config_path);
    json program_args = json::parse(f);

    for (auto& experiment : program_args["experiments"]) {
      auto simulation = Simulation(experiment);
      std::cout << "Starting experiment " << experiment["name"] << std::endl;
      if (experiment.contains("opath")) {
        TIMERSTART(sim)
        simulation.start(experiment["tot_time"],
                         experiment["time_step"],
                         experiment["opath"]);
        TIMERSTOP(sim)
      } else {
        TIMERSTART(sim)
        simulation.start(experiment["tot_time"], experiment["time_step"]);
        TIMERSTOP(sim)
      }
    }
  }

  // -xAVX2 and -xMIC-AVX512 flags force the compiler to generate AVX2
  // and AVX-512 SIMD instructions, respectively. AVX2 extensions accelerated
  // the previous version by a factor of 7.4× while AVX-512 instructions
  // achieved a speedup of 15.1×
  // http://sedici.unlp.edu.ar/bitstream/handle/10915/95855/Documento_completo.pdf?sequence=1

  return 0;
}