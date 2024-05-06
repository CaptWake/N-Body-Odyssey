# Define directories
SIM_HEADER_DIR = include/simulations/
UTILS_DIR = include/utilities/
SRC_DIR = src/simulations/

# Define C++ compiler and flags
CXX = g++
CXXFLAGS = -Wall -O3 -std=c++17
AVX_FLAGS = -march=native
OPENMP_FLAGS = -fopenmp

# Define object files pattern
OBJECTS = $(SRC_DIR)%.o

# Define all available executables
ALL_EXECUTABLES = nbody_sequential_ap nbody_sequential_ap_avx nbody_omp_ap nbody_sequential_bh nbody_sequential_bh_avx nbody_omp_bh

.PHONY: all clean

all: $(ALL_EXECUTABLES)

nbody_sequential_ap: $(SRC_DIR)sequential_ap.cc
	$(CXX) $(CXXFLAGS) $(SRC_DIR)sequential_ap.cc -I$(SIM_HEADER_DIR) -I$(UTILS_DIR) -o $@

nbody_sequential_ap_avx: $(SRC_DIR)sequential_ap_avx.cc
	$(CXX) $(CXXFLAGS) $(AVX_FLAGS) $(SRC_DIR)sequential_ap_avx.cc -I$(SIM_HEADER_DIR) -I$(UTILS_DIR) -o $@

nbody_omp_ap: $(SRC_DIR)omp_ap.cc
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(SRC_DIR)omp_ap.cc -I$(SIM_HEADER_DIR) -I$(UTILS_DIR) -o $@ -DOMP

nbody_sequential_bh: $(SRC_DIR)sequential_bh.cc
	$(CXX) $(CXXFLAGS) $(SRC_DIR)sequential_bh.cc -I$(SIM_HEADER_DIR) -I$(UTILS_DIR) -o $@

nbody_sequential_bh_avx: $(SRC_DIR)sequential_bh_avx.cc
	$(CXX) $(CXXFLAGS) $(AVX_FLAGS) $(SRC_DIR)sequential_bh_avx.cc -I$(SIM_HEADER_DIR) -I$(UTILS_DIR) -o $@

nbody_omp_bh: $(SRC_DIR)omp_bh.cc
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(SRC_DIR)omp_bh.cc -I$(SIM_HEADER_DIR) -I$(UTILS_DIR) -o $@ -DOMP

# Generic object file compilation rule
$(OBJECTS): $(SRC_DIR)%.cc
	$(CXX) $(CXXFLAGS) $< -I$(SIM_HEADER_DIR) -I$(UTILS_DIR) -o $@

clean:
	rm -f $(ALL_EXECUTABLES) $(OBJECTS)

