#include <argparse/argparse.hpp>
#include <fstream>

#include "simulation.h"
#include "time_utils.h"
#include "octree.h"

int main(int argc, char** argv) {
  std::vector<float> px {0.2f, 0.1f, 0.3f, 0.8f, 0.971f};
  std::vector<float> py {0.71f, 0.26f, 0.14f, 0.55f, 0.34f};
  std::vector<float> pz {0.28f, 0.444f, 0.645f, 0.222f, 0.04f};
  std::vector<float> m {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  auto p = octreeSOA(px, py, pz, m);
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