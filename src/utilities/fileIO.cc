#include "fileIO.h"

#include <fstream>
#include <iostream>
#include <sstream>

uint64_t ReadCSVConfigurationAOS(const std::string& filename,
                                 std::vector<float>& m, std::vector<float>& p,
                                 std::vector<float>& v, float& G) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
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

uint64_t ReadCSVConfigurationSOA(const std::string& filename,
                                 std::vector<float>& m, std::vector<float>& px,
                                 std::vector<float>& py, std::vector<float>& pz,
                                 std::vector<float>& vx, std::vector<float>& vy,
                                 std::vector<float>& vz, float& G) {
  std::ifstream file(filename);
  if (file.is_open()) {
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
      float mass, x, y, z, vx_, vy_, vz_;
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
      vx_ = std::stof(token);

      std::getline(iss, token, ',');
      vy_ = std::stof(token);

      std::getline(iss, token, ',');
      vz_ = std::stof(token);

      m.push_back(mass);
      px.push_back(x);
      py.push_back(y);
      pz.push_back(z);
      vx.push_back(vx_);
      vy.push_back(vy_);
      vz.push_back(vz_);
      i++;
    }
    file.close();
    return i;
  } else {
    std::cerr << "Unable to open file: " << filename << std::endl;
    exit(1);
  }
}