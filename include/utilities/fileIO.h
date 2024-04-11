#ifndef FILEIO_H
#define FILEIO_H

#include <vector>
#include <cstdint>
#include <string>
#include <immintrin.h>

uint64_t ReadCSVConfigurationAOS(const std::string &filename, std::vector<float>& m, std::vector<float>& p, std::vector<float>& v, float& G);
uint64_t ReadCSVConfigurationSOA(const std::string &filename, std::vector<float>& m, std::vector<float>& px, std::vector<float>& py, std::vector<float>& pz, std::vector<float>& vx, std::vector<float>& vy, std::vector<float>& vz, float& G);

template <typename T>
void vector_to_arr(std::vector<T> vec, T **arr, bool use_avx) {
  if (use_avx) {
    // rounds the number of slots to multiple of 8
    uint64_t n_slots = 8 * floor(vec.size() / 8 + 1);
    *arr = static_cast<T*>(_mm_malloc(n_slots * sizeof(T), 32));
  }
  else
    *arr = new float[vec.size()];
   std::copy(vec.begin(), vec.end(), *arr);
}

#endif