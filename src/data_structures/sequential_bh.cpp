#include "sequential_bh.h"
#include "octree.h"
#include <fstream>
#include <sstream>

void SequentialBH::Update(float dt) {
  std::vector<vec3> f{this->p.size(), {0,0,0}};
  octree tree(this->p, this->m);
  for (auto i = 0; i < this->p.size(); ++i){
    f[i] = tree.force_at(this->p[i], 0, theta);
    this->v[i] += f[i] * dt;
  }

  for (std::size_t i = 0; i < p.size(); ++i) {
    this->p[i] += this->v[i] * dt;
  }
}

void SequentialBH::LoadFromCSVConfiguration(const std::string &filename) {
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
  this->G = std::stof(token);

  std::getline(file, line);
  iss = std::istringstream(line);
  std::getline(iss, token);
  // Get the theta parameter (second line)
  this->theta = std::stof(token);


  float mass, x, y, z, vx, vy, vz;
  while (std::getline(file, line)) {
    iss = std::istringstream(line);
    // Read particle data from CSV fields parsing comma-separated values
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

    this->m.push_back(mass);
    this->p.emplace_back(x,y,z);
    this->v.emplace_back(vx, vy, vz);
    i++;
  }
  file.close();
}

void SequentialBH::LogsToCSV(const std::string &filename) {
  std::ofstream file(filename, std::ios_base::app);
  if (file.is_open()) {
    auto n_bodies = this->p.size();
    for (uint64_t i = 0; i < n_bodies; ++i) {
      file << this->p[i].x() << "," << this->p[i].y() << "," << this->p[i].z();
      if (i < n_bodies - 1)
        file << ",";
    }
    file << std::endl;
    file.close();
  } else {
    std::cerr << "Unable to open file: " << filename << std::endl;
  }
}

