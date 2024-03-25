#include "../include/simulation.h"

int main() {
  auto simulation = Simulation();
  simulation.start(100, 0.01);
  return 0;
}