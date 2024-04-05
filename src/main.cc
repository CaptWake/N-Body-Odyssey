#include "simulation.h"
#include "time_utils.h"
#include <argparse/argparse.hpp>

#include "quadtree.h"

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
  //simulation = Simulation("SEQ", R"(D:\universita\SCPD\project1\input.csv)");
  //TIMERSTART(SEQUENTIAL)
  //simulation.start(10, 0.01, "");
  //TIMERSTOP(SEQUENTIAL)

  // -xAVX2 and -xMIC-AVX512 flags force the compiler to generate AVX2
  //and AVX-512 SIMD instructions, respectively. AVX2 extensions accelerated the previous version by a factor of 7.4×
  //while AVX-512 instructions achieved a speedup of 15.1×
  // http://sedici.unlp.edu.ar/bitstream/handle/10915/95855/Documento_completo.pdf?sequence=1
  //simulation = Simulation("SEQ_AVX", program.get<std::string>("-i"));
  //TIMERSTART(SEQUENTIAL_AVX)
  //simulation.start(10, 0.01, "/home/ste/Documents/test2.csv");
  //TIMERSTOP(SEQUENTIAL_AVX)
  std::vector<float> m;
  std::vector<vec3> p;
  std::vector<vec3> v;
  LoadFromCSVConfiguration(R"(D:\universita\SCPD\project1\input.csv)", m, p, v);
  for (int i = 0; i < 300; ++i) {
    update(p, m, v, 0.01f);
    std::cout << p[0] << "; " << p[1] << "; " << p[2] << std::endl;
  }

  return 0;
}