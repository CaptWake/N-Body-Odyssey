#ifndef DATA_FORMAT_H_
#define DATA_FORMAT_H_

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

// Function to export bodies information to CSV file
void ExportToCSV(const float* data, const uint64_t n_bodies,
                 const std::string& filename) {
  std::ofstream file(filename, std::ios_base::app);
  if (file.is_open()) {
    uint64_t i;
    for (i = 0; i < n_bodies*3; i+=3) {
      std::cout << data[i] << "," << data[i+1] << "," << data[i+2] << std::endl;
      file << data[i] << "," << data[i+1] << "," << data[i+2];
      if (i < n_bodies * 3 - 3)
        file << ",";
    }
    file << std::endl;
    file.close();
  } else {
    std::cerr << "Unable to open file: " << filename << std::endl;
  }
}

uint64_t LoadFromCSVConfiguration(const std::string& filename, float** masses, float** positions, float** velocities, float& grav_const) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    return 0;  // Indicate error by returning 0 particles loaded
  }

  std::string line;
  uint64_t i = 0;

  // Get header line (optional, but recommended for clarity)
  // std::getline(file, line);  // Assuming header is present

  std::vector<float> m, p, v;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string token;

    // Get the gravitational constant (first line)
    if (i == 0) {
      std::getline(iss, token);
      grav_const = std::stof(token);
    }
    else {
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
    }
    i++;
  }

  i--;

  file.close();

  // Update dynamic arrays
  *masses = (float *)malloc(i * sizeof(float));
  *positions = (float *)malloc(3*i * sizeof(float));
  *velocities = (float *) malloc(3*i * sizeof(float));

  for (uint64_t k = 0; k < 3*i; ++k) {
    std::cout << "Mass: " << m[k/3] << std::endl;
    (*masses)[k/3] = m[k/3];
    std::cout << "Position: " << p[k] << std::endl;
    (*positions)[k] = p[k];
    std::cout << "Velocity: " << v[k] << std::endl;
    (*velocities)[k] = v[k];
  }
  std::cout << "GRAVITY: " << grav_const << std::endl;
  std::cout << "NBODIES: " << i << std::endl;
  return i;  // Return the number of particles loaded
}


#endif
