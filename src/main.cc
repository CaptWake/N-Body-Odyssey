#include "simulation.h"
#include "time_utils.h"
#include <argparse/argparse.hpp>

int main(int argc, char** argv) {
  argparse::ArgumentParser program("N-Body-Simulator");
  program.add_argument("-a", "--algorithm")
      .default_value(std::string{"AP"})
      .choices("AP", "BH")
      .help("specify the algorithm");
  program.add_argument("-m", "--mode")
      .default_value(std::string{"SEQ"})
      .choices("SEQ", "SEQ_AVX", "CUDA", "OMP")
      .help("specify the simulation mode");
  program.add_argument("-i", "--input")
      .default_value(std::string{""})
      .help("specify the input data file path");
  program.add_argument("-o", "--output")
      .default_value(std::string{""})
      .help("specify the output file path where simulation results will be stored");

  try {
    program.parse_args(argc, argv);    // Example: ./main --color orange
  }
  catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  Simulation simulation;

  //if (program.is_used("-i")) {
  //  simulation = Simulation(program.get<std::string>("-m"), program.get<std::string>("-i"));
  //}
  simulation = Simulation("SEQ", program.get<std::string>("-i"));
  TIMERSTART(SEQUENTIAL)
  simulation.start(10, 0.01);
  TIMERSTOP(SEQUENTIAL)

  // -xAVX2 and -xMIC-AVX512 flags force the compiler to generate AVX2
  //and AVX-512 SIMD instructions, respectively. AVX2 extensions accelerated the previous version by a factor of 7.4×
  //while AVX-512 instructions achieved a speedup of 15.1×
  // http://sedici.unlp.edu.ar/bitstream/handle/10915/95855/Documento_completo.pdf?sequence=1
  simulation = Simulation("SEQ_AVX", program.get<std::string>("-i"));
  TIMERSTART(SEQUENTIAL_AVX)
  simulation.start(10, 0.01);
  TIMERSTOP(SEQUENTIAL_AVX)

  return 0;
}