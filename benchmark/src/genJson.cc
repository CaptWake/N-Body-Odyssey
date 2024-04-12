#include <iostream>
#include <string>

#define N_ALGOS 7
#define N_TEST 15

int main(int argc, char** argv) {
  std::string name[] = {"All Pairs Sequential",
                        "All Pairs Sequential with AVX instructions",
                        "All Pairs with OpenMP - static schedule",
                        "All Pairs with OpenMP - dynamic schedule",
                        "Barnes-Hut Sequential",
                        "Barnes-Hut Sequential - dynamic schedule",
                        "Barnes-Hut Sequential - static schedule"};
  std::string algo[] = {"AP", "AP", "AP", "AP", "BH", "BH", "BH"};
  std::string mode[] = {"SEQ", "SEQ_AVX", "OMP", "OMP", "SEQ", "OMP", "OMP"};
  float step = 0.01;
  float tot_time = 10;
  int n_bodies = 2;
  float G = 1;
  std::cout << "{\n";
  std::cout << "\t\"experiments\": [\n";
  for (int j = 0; j < N_TEST; j++) {
    for (int i = 0; i < N_ALGOS; i++) {
      std::cout << "\t\t{\n";
      std::cout << "\t\t\t\"name\": \"" << name[i] << "\",\n";
      std::cout << "\t\t\t\"algorithm\": \"" << algo[i] << "\",\n";
      std::cout << "\t\t\t\"mode\": \"" << mode[i] << "\",\n";
      std::cout << "\t\t\t\"time_step\": " << step << ",\n";
      std::cout << "\t\t\t\"tot_time\": " << tot_time << ",\n";
      std::cout << "\t\t\t\"n_bodies\": " << n_bodies << ",\n";
      if (i == 2 || i == 5 || i == 3 || i == 6) {
        std::cout << "\t\t\t\"num_threads\": 8,\n";
        std::cout << "\t\t\t\"chunk_size\": 1,\n";
      }
      if (i == 2 || i == 5) {
        std::cout << "\t\t\t\"schedule_type\": \"static\",\n";
      }
      if (i == 3 || i == 6) {
        std::cout << "\t\t\t\"schedule_type\": \"dynamic\",\n";
      }
      if (i == 4) {
        std::cout << "\t\t\t\"theta\": 0.9,\n";
      }
      std::cout << "\t\t\t\"grav_const\": " << G << "\n";
      std::cout << "\t\t},\n";
    }
    n_bodies *= 2;
  }
  std::cout << "\t]\n";
  std::cout << "}\n";
  return 0;
}