#include "../include/simulation.h"

int main() {
  auto simulation = SequentialAllPairsSimulation();
  simulation.start(10, 0.1);
  return 0;
}