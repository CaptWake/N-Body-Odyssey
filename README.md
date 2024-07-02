# N-Body Odyssey: Algorithms and Technologies

This project explores N-Body simulations using different algorithms and technologies to study their performance characteristics. It includes implementations of All-Pairs and Barnes-Hut algorithms, leveraging various technologies such as OpenMP, OpenMPI, Rust, CUDA, and AVX.

## Project Structure

The project structure is organized as follows:

- `format.sh`: Shell script for formatting source code.
- `include/`: Directory for header files.
  - `simulations/`: Header files for simulation algorithms.
  - `utilities/`: Utility header files.
- `Makefile`: Makefile for building the project.
- `src/`: Source code directory.
  - `simulations/`: Source files for simulation algorithms.
  - `utilities/`: Source files for utility functions.
- `tests/`: Directory for test-related files.

## Technologies Used

The project utilizes several technologies to explore different aspects of N-Body simulation:

- **OpenMP**: Used for shared-memory parallelization on multicore CPUs.
- **OpenMPI**: Employed for distributed-memory parallelization across multiple nodes or processors.
- **Rust**: Chosen for its popularity about performance, memory safety, and modern language features.
- **CUDA**: NVIDIA's parallel computing platform and programming model used for GPU acceleration.
- **AVX**: Advanced Vector Extensions used to enhance performance via SIMD (Single Instruction, Multiple Data) operations on compatible intel CPUs.

## Algorithms Implemented

### All-Pairs Algorithm `O(n^2)`

- **Sequential Implementation**: Basic implementation without parallelism.
- **OpenMP**: Parallel version using OpenMP directives to exploit multicore processors.
- **OpenMPI**: Implementation using MPI (Message Passing Interface) for distributed memory parallelism across multiple nodes.
- **AVX**: Implementation leveraging AVX (Advanced Vector Extensions) for SIMD (Single Instruction, Multiple Data) acceleration.
- **OpenMPI + OpenMP + AVX**: Hybrid implementation combining MPI for distributed computing, OpenMP for shared-memory parallelism within nodes, and AVX for SIMD acceleration.

### Barnes-Hut Algorithm `O(n log n)`

- **Sequential Implementation**: Hierarchical algorithm for reducing computational complexity by approximating distant groups of particles as single entities.

## How to Run Experiments

To run the experiments, follow these steps:

1. **Build the Project**:
   - Navigate to the project directory containing the `Makefile`.
   - Use the `make` command to build the project. This will compile all executable files.

2. **Run Specific Experiments**:
   - Choose the executable based on the desired algorithm and technology combination.
   - Execute the chosen executable with appropriate command-line arguments (if any).

Example:
```bash
make nbody_sequential_ap_avx 
./nbody_sequential_ap_avx <number_of_particles>
```

Replace `<number_of_particles>` with the desired number of particles for the simulation.

3. **Evaluate Results**:
   - Monitor console output for simulation results, including performance metrics or any errors encountered.

## Notes

- Adjust configurations or parameters in `Makefile` to tailor simulations or experiments.
- Ensure dependencies for each technology (e.g., CUDA drivers, MPI libraries) are properly installed and configured.

