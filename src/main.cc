#include "../include/simulation.h"

int main() {
  auto simulation = Simulation();
  simulation.start(0.03, 0.01);
  return 0;
}