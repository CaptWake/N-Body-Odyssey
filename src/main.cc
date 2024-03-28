#include "../include/simulation.h"
#include <argparse/argparse.hpp>

int main(int argc, char** argv) {
  argparse::ArgumentParser program("test");
  program.add_argument("-a", "--algorithm")
      .default_value(std::string{"AP"})
      .choices("AP", "BH")
      .help("specify the algorithm");
  program.add_argument("-m", "--mode")
      .default_value(std::string{"SEQ"})
      .choices("SEQ", "CUDA", "OMP")
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

  if (program.is_used("-i")) {
    simulation = Simulation(program.get<std::string>("-i"));
  }
  simulation.start(100, 0.01, program.get<std::string>("-o"));
  //auto color = program.get<std::string>("--color");  // "orange"
  //auto explicit_color = program.is_used("--color");  // true, user provided orange
  //std::cout << color << std::endl;

  //auto simulation = Simulation();
  //simulation.start(100, 0.01);
  return 0;
}