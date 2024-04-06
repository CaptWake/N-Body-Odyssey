#include "simulation.h"
#include "time_utils.h"
#include <argparse/argparse.hpp>

#include "octree.h"

// Recursive function to print each node
void printNodes(const node& currentNode, const node_id& currentId)
{
  // Print the current node ID
  std::cout << "Node ID: " << currentId << std::endl;

  // Iterate over children and recursively print them
  for (int i = 0; i < 8; ++i) {
        node_id childId = currentNode.children[i];
        if (childId != null) {
          // Print the child node
          std::cout << "Child[" << i << "]: " << childId << std::endl;
        }
    }
}

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
  simulation = Simulation("AP", "SEQ", "");
  TIMERSTART(SEQUENTIAL_AP)
  simulation.start(10, 0.01);
  TIMERSTOP(SEQUENTIAL_AP)

  // -xAVX2 and -xMIC-AVX512 flags force the compiler to generate AVX2
  //and AVX-512 SIMD instructions, respectively. AVX2 extensions accelerated the previous version by a factor of 7.4×
  //while AVX-512 instructions achieved a speedup of 15.1×
  // http://sedici.unlp.edu.ar/bitstream/handle/10915/95855/Documento_completo.pdf?sequence=1
  simulation = Simulation("AP", "SEQ_AVX", "");
  TIMERSTART(SEQUENTIAL_AP_AVX)
  simulation.start(10, 0.01);
  TIMERSTOP(SEQUENTIAL_AP_AVX)

  simulation = Simulation("BH", "SEQ", "");
  TIMERSTART(SEQUENTIAL_BH)
  simulation.start(10, 0.01);
  TIMERSTOP(SEQUENTIAL_BH)

  return 0;
}