# Define directories

INCLUDE_DIR = -I include/
SRC_DIR = src/simulations/

# Define C++ compiler and flags
CXX = g++
MPIPPC = mpic++
NVCXX = nvcc 
CXXFLAGS = -Wall -O3 -std=c++17
AVX_FLAGS = -march=native
OPENMP_FLAGS = -fopenmp
PRECISION = FLOAT
# Define object files pattern
OBJECTS = $(SRC_DIR)%.o

# GPU settings
GPU_RUNTIME := CUDA

# HIP variables
ROCM_INSTALL_DIR := /opt/rocm
HIP_INCLUDE_DIR  := $(ROCM_INSTALL_DIR)/include

HIPCXX ?= $(ROCM_INSTALL_DIR)/bin/hipcc

# Common variables and flags
CXX_STD   := c++17
ICXXFLAGS := -std=$(CXX_STD)
ICPPFLAGS := $(COMMON_INCLUDE_DIR)
ILDFLAGS  :=
ILDLIBS   :=

ifeq ($(GPU_RUNTIME), CUDA)
	ICXXFLAGS += -x cu
	ICPPFLAGS += -isystem $(HIP_INCLUDE_DIR)
else ifeq ($(GPU_RUNTIME), HIP)
	CXXFLAGS ?= -Wall -Wextra
else
	$(error GPU_RUNTIME is set to "$(GPU_RUNTIME)". GPU_RUNTIME must be either CUDA or HIP)
endif

ICXXFLAGS += $(CXXFLAGS)
ICPPFLAGS += $(CPPFLAGS)
ILDLIBS   += $(LDLIBS)

# Define all available executables
ALL_EXECUTABLES = nbody_sequential_ap nbody_sequential_ap_avx nbody_omp_ap nbody_sequential_bh nbody_sequential_bh_avx nbody_omp_bh nbody_hip_ap nbody_cuda_ap_double nbody_cuda_ap_float nbody_mpi_ap nbody_mpi_omp_avx_ap

.PHONY: all clean

all: $(ALL_EXECUTABLES)

nbody_sequential_ap: $(SRC_DIR)sequential_ap.cc
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(SRC_DIR)sequential_ap.cc $(INCLUDE_DIR) -o $@ -D$(PRECISION) -DOMP #-DMONITOR_ENERGY 

nbody_sequential_ap_avx: $(SRC_DIR)sequential_ap_avx.cc
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(AVX_FLAGS) $(SRC_DIR)sequential_ap_avx.cc $(INCLUDE_DIR) -o $@ -D$(PRECISION) -DOMP -DMONITOR_ENERGY 

nbody_omp_ap: $(SRC_DIR)omp_ap.cc
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(SRC_DIR)omp_ap.cc $(INCLUDE_DIR) -o $@ -DOMP -D$(PRECISION) #-DMONITOR_ENERGY

nbody_mpi_ap: $(SRC_DIR)mpi_ap.cc
	$(MPIPPC) $(CXXFLAGS) $(OPENMP_FLAGS) $(SRC_DIR)mpi_ap.cc $(INCLUDE_DIR) -o $@ -D$(PRECISION) -DOMP -DUSE_MPI -DMONITOR_ENERGY

nbody_sequential_bh: $(SRC_DIR)sequential_bh.cc
	$(CXX) $(CXXFLAGS) $(SRC_DIR)sequential_bh.cc $(INCLUDE_DIR) -o $@

nbody_sequential_bh_avx: $(SRC_DIR)sequential_bh_avx.cc
	$(CXX) $(CXXFLAGS) $(AVX_FLAGS) $(SRC_DIR)sequential_bh_avx.cc $(INCLUDE_DIR) -o $@

nbody_omp_bh: $(SRC_DIR)omp_bh.cc
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(SRC_DIR)omp_bh.cc $(INCLUDE_DIR) -o $@ -DOMP

nbody_mpi_omp_avx_ap: $(SRC_DIR)mpi_omp_avx_ap.cc
	$(MPIPPC) $(CXXFLAGS) $(OPENMP_FLAGS) $(AVX_FLAGS) -DOMP $(SRC_DIR)mpi_omp_avx_ap.cc $(INCLUDE_DIR) -o $@ -D$(PRECISION) -DOMP -DUSE_MPI -DMONITOR_ENERGY

nbody_hip_ap: $(SRC_DIR)ap_soa.hip
	$(HIPCXX) $(ICXXFLAGS) $(ICPPFLAGS) $(ILDFLAGS) $(INCLUDE_DIR) -o $@ $< $(ILDLIBS) -D$(PRECISION)

nbody_cuda_ap_float: $(SRC_DIR)ap_soa_float.cu
	$(NVCXX) $(INCLUDE_DIR) -Xcompiler "-fopenmp" -o $@ $< 

nbody_cuda_ap_double: $(SRC_DIR)ap_soa_double.cu
	$(NVCXX) $(INCLUDE_DIR) -Xcompiler "-fopenmp" -o $@ $< 

# Generic object file compilation rule
$(OBJECTS): $(SRC_DIR)%.cc
	$(CXX) $(CXXFLAGS) $< $(INCLUDE_DIR) -o $@

clean:
	rm -f $(ALL_EXECUTABLES) $(OBJECTS)
