#ifndef DATA_FORMAT_H_
#define DATA_FORMAT_H_

#include <fstream>
#include <iostream>
#include <string>

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

#endif
