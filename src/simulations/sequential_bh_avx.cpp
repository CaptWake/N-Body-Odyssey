#include "sequential_bh_avx.h"
#include "octree.h"
#include "avx.h"
#include <fstream>
#include <sstream>


void SequentialBHAVX::Update(float dt) {
  node_id children[8]{0, null, null, null, null, null, null, null};

  auto tree = octreeSOA(this->px, this->py, this->pz, this->m, this->n_bodies);
  for (uint64_t i = 0; i < this->n_bodies; ++i) {
    __m256 Fx = _mm256_set1_ps(0.f);
    __m256 Fy = _mm256_set1_ps(0.f);
    __m256 Fz = _mm256_set1_ps(0.f);

    tree.force_at(this->px[i], this->py[i], this->pz[i], children, theta, Fx, Fy, Fz);

    this->vx[i] += hsum_avx(Fx) * dt;
    this->vy[i] += hsum_avx(Fy) * dt;
    this->vz[i] += hsum_avx(Fz) * dt;

  }

  __m256 Dt = _mm256_set1_ps(dt);

  for (uint64_t i = 0; i < this->n_bodies; i+=8) {
    __m256 Px = _mm256_load_ps(this->px+i);
    __m256 Py = _mm256_load_ps(this->py+i);
    __m256 Pz = _mm256_load_ps(this->pz+i);

    __m256 Vx = _mm256_load_ps(this->vx+i);
    __m256 Vy = _mm256_load_ps(this->vy+i);
    __m256 Vz = _mm256_load_ps(this->vz+i);

    _mm256_store_ps(this->px+i, _mm256_fmadd_ps(Vx, Dt, Px));
    _mm256_store_ps(this->py+i, _mm256_fmadd_ps(Vy, Dt, Py));
    _mm256_store_ps(this->pz+i, _mm256_fmadd_ps(Vz, Dt, Pz));

  }
}

void SequentialBHAVX::LogsToCSV(const std::string &filename) {
  std::ofstream file(filename, std::ios_base::app);
  if (file.is_open()) {
    for (uint64_t i = 0; i < this->n_bodies; ++i) {
      file << this->px[i] << "," << this->py[i] << "," << this->pz[i];
      if (i < n_bodies - 1) file << ",";
    }
    file << std::endl;
    file.close();
  } else {
    std::cerr << "Unable to open file: " << filename << std::endl;
  }
}
