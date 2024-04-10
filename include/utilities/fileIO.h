#ifndef FILEIO_H
#define FILEIO_H

#include <vector>
#include <cstdint>
#include <string>

uint64_t ReadCSVConfiguration(const std::string &filename, std::vector<float>& m, std::vector<float>& p, std::vector<float>& v, float& G);

#endif