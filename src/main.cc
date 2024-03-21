#include "../include/simulation.h"

int main() {
  auto simulation = Simulation();
  simulation.start(10, 0.1);
  return 0;
}