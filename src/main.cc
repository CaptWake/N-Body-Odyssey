#include "../data_structures/particle.h"

int main() {
  auto mercury = Particle(1.6601367952719304e7,
                          vec3(3.067925455119497e1, -2.695851754307215e-01, -5.017408038016835e-02),
                          vec3(1.300791785131080e-02, 2.245507195129482e-02, 6.414809034847098e-04),
                          {0, 0, 0}, 20, 5.43);

  std::cout << mercury << std::endl;
  return 0;
}