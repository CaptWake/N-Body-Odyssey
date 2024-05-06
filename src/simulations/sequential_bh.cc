#include "sequential_bh.h"

#include <fstream>
#include <sstream>

#include "octree.h"

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

void SequentialBH::LogsToCSV(const std::string &filename) {
  std::ofstream file(filename, std::ios_base::app);
  if (file.is_open()) {
    auto n_bodies = this->p.size();
    for (uint64_t i = 0; i < n_bodies; ++i) {
      file << this->p[i].x() << "," << this->p[i].y() << "," << this->p[i].z();
      if (i < n_bodies - 1) file << ",";
    }
    file << std::endl;
    file.close();
  } else {
    std::cerr << "Unable to open file: " << filename << std::endl;
  }
}
