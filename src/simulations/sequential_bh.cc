#include "simulations/sequential_bh.h"
#include "utilities/octree.h"
#include "utilities/time_utils.h"

void SequentialBH::Update(float dt) {
  std::vector<vec3> f{this->p.size(), {0, 0, 0}};
  octree tree(this->p, this->m);
  for (auto i = 0; i < this->p.size(); ++i) {
    f[i] = tree.force_at(this->p[i], 0, theta);
    this->v[i] += f[i] * dt;
  }

  for (std::size_t i = 0; i < p.size(); ++i) {
    this->p[i] += this->v[i] * dt;
  }
}

float Ep(uint64_t n, const float *m, const vec3 *p) {
  float Epot = 0.0;
  uint64_t i, j;
  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
      auto d = p[i] - p[j];
      auto D = d.length();
      Epot += -1.0 * m[i] * m[j] / D;
    }
  }
  return Epot;
}

float Ek(uint64_t n, const float *m, const vec3 *v) {
  float Ekin = 0.0;
  uint64_t i;
  for (i = 0; i < n; ++i) {
    Ekin += 0.5 * m[i] * v[i].squared_length();
  }
  return Ekin;
}

void SequentialBHSimulate(uint64_t n, float theta, float dt, float tEnd) {
  auto simulation = SequentialBH(n, theta);
  for (float t = 0; t < tEnd; t+=dt) {
    simulation.Update(dt);
    auto ek = Ek(n, simulation.m.data(), simulation.v.data());
    auto ep = Ep(n, simulation.m.data(), simulation.p.data());
    std::cout << " Etot: " << ek + ep << std::endl;
  }
}
//
int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Must specify the number of bodies and theta" << std::endl;
    exit(1);
  }
  srand(0);
  TIMERSTART(simulation)
  SequentialBHSimulate(std::stoul(argv[1]), 0.1, 0.01, 0.1);
  TIMERSTOP(simulation)
}