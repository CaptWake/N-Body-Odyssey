#include "fileIO.h"

#include <fstream>
#include <sstream>
#include <iostream>

uint64_t ReadCSVConfiguration(const std::string &filename, std::vector<float>& m, std::vector<float>& p, std::vector<float>& v, float& G) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Cannot find configuration file at: " << filename << std::endl;
    exit(1);
  }

  std::string line;
  uint64_t i = 0;

  // Assuming the header is not present

  std::getline(file, line);
  std::istringstream iss(line);
  std::string token;
  std::getline(iss, token);
  // Get the gravitational constant (first line)
  G = std::stof(token);

  while (std::getline(file, line)) {
    iss = std::istringstream(line);
    // Read particle data from CSV fields
    float mass, x, y, z, vx, vy, vz;
    // Read and parse comma-separated values
    std::getline(iss, token, ',');
    mass = std::stof(token);

    std::getline(iss, token, ',');
    x = std::stof(token);

    std::getline(iss, token, ',');
    y = std::stof(token);

    std::getline(iss, token, ',');
    z = std::stof(token);

    std::getline(iss, token, ',');
    vx = std::stof(token);

    std::getline(iss, token, ',');
    vy = std::stof(token);

    std::getline(iss, token, ',');
    vz = std::stof(token);

    m.push_back(mass);
    p.push_back(x);
    p.push_back(y);
    p.push_back(z);
    v.push_back(vx);
    v.push_back(vy);
    v.push_back(vz);
    i++;
  }
  file.close();

  return i;
}