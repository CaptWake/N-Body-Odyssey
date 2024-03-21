#ifndef DATA_FORMAT_H_
#define DATA_FORMAT_H_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <particle.h>
#include "vec3.h"

// Function to export particles information to CSV file
void ExportToCSV(const std::vector<Particle>& data,
                 const std::string& filename) {
  std::ofstream file(filename, std::ios_base::app);
  if (file.is_open()) {
    auto nParticles = data.size();
    for (int i = 0; i < nParticles; ++i) {
      auto point = data[i].GetPosition();
      file << point.x() << "," << point.y() << "," << point.z();
      if (i == nParticles - 1)
        file << std::endl;
      else
        file << ",";
    }
    file.close();
  } else {
    std::cerr << "Unable to open file: " << filename << std::endl;
  }
}

#endif
