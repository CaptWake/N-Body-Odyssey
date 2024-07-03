#include "utilities/time_utils.h"
#include "utilities/nbody_helpers.h"
#include <iostream>
#include <omp.h>

#define MAX_THREADS_PER_BLOCK 1024
#include <cstring> // Add this line to include the necessary header for strcmp
#include <cstdlib> // Add this line to include the necessary header for exit
#include <algorithm>

/**
 * Sets the CUDA device based on the provided device name.
 *
 * @param deviceName The name of the device to set.
 */
void setDeviceByName(const char* deviceName) {
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  int targetDevice = -1;
  for (int i = 0; i < deviceCount; ++i) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, i);
    std::string fullName(deviceProp.name);
    std::string shortName(deviceName);
    // Convert both names to lowercase for case-insensitive comparison
    std::transform(fullName.begin(), fullName.end(), fullName.begin(), ::tolower);
    std::transform(shortName.begin(), shortName.end(), shortName.begin(), ::tolower);
    if (fullName.find(shortName) != std::string::npos) {
      targetDevice = i;
      break;
    }
  }
  if (targetDevice != -1) {
    cudaSetDevice(targetDevice);
    std::cout << "Device " << deviceName << " set as the current device." << std::endl;
  } else {
    std::cerr << "Device " << deviceName << " not found." << std::endl;
    exit(1);
  }
}

double inline xoxyi_rand(unsigned int *seed){
  return (double)rand_r(seed) / (double) RAND_MAX;
}

/**
 * Calculates the potential energy (Epot) of a system of particles using the Lennard-Jones potential.
 *
 * @param n The number of particles in the system.
 * @param p Pointer to an array of double4 structures representing the particles.
 *          Each double4 structure contains the x, y, z coordinates of a particle (p[i].x, p[i].y, p[i].z)
 *          and the weight of the particle (p[i].w).
 * @return The total potential energy of the system.
 */
double inline Ep(const int n, double4 *p) {
  double Epot = 0.0f;
#pragma omp parallel for reduction(+:Epot)
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double D, x, y, z;
      x = p[i].x - p[j].x;
      y = p[i].y - p[j].y;
      z = p[i].z - p[j].z;
      D = sqrt(x * x + y * y + z * z);
      Epot += -1.0f * p[i].w * p[j].w / D;
    }
  }
  return Epot;
}

/**
 * Calculates the kinetic energy of a system of particles.
 *
 * @param n The number of particles in the system.
 * @param v Pointer to an array of double4 structures representing the particles' properties.
 *           Each double4 structure contains the particle's position (x, y, z) and mass (w).
 * @return The total kinetic energy of the system.
 */
double inline Ek(const int n, double4 *v) {
  double Ekin = 0.0;
#pragma omp parallel for reduction(+:Ekin)
  for (int i = 0; i < n; ++i) {
    Ekin += 0.5f * v[i].w * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z);
  }
  return Ekin;
}

/**
 * Scales a 4D array by a given scale factor.
 *
 * @param n     The size of the array.
 * @param m     The 4D array to be scaled.
 * @param scale The scale factor to multiply each element of the array by.
 */
static inline void scale4DArray(const int n, double4 *m, const double scale) {
#pragma omp parallel for 
  for (int i = 0; i < n; ++i) {
    m[i].x *= scale;
    m[i].y *= scale;
    m[i].z *= scale;
  }
}

/**
 * @brief Initializes the position and velocity arrays for a given number of particles.
 * 
 * @param n The number of particles.
 * @param p Pointer to the position array.
 * @param v Pointer to the velocity array.
 * @param seed The seed value for the random number generator.
 */
void InitPosVel(const int n, double4 *p, double4 *v, int seed) {
#pragma omp parallel
  {
    unsigned int seedT = omp_get_thread_num() + seed * omp_get_num_threads();	 
    #pragma omp for
    for (int i = 0; i < n; ++i) {
      double R, X, Y;
      double mi = _M / n;
      R = xoxyi_rand(&seedT);
      X = acos(1.0f - 2.0f * xoxyi_rand(&seedT));
      Y = xoxyi_rand(&seedT) * 2.0f * _PI;
      // https://www.researchgate.net/figure/Figure-A1-Spherical-coordinates_fig8_284609648
      p[i] = make_double4(R * sin(X) * cosf(Y), R * sinf(X) * sinf(Y), R * cosf(X), mi);
      v[i] = make_double4(1.0f - 2.0f * xoxyi_rand(&seedT), 1.0f - 2.0f * xoxyi_rand(&seedT), 1.0f - 2.0f * xoxyi_rand(&seedT), mi);
    }
  }
}

/**
 * Moves the particles to the center of mass and subtracts the center of mass velocity from each particle.
 *
 * @param n The number of particles.
 * @param p Pointer to an array of double4 structures representing the position of each particle.
 * @param v Pointer to an array of double4 structures representing the velocity of each particle.
 */
void Move2Center(const int n, double4 *p, double4 *v) {
  double3 pp = make_double3(0.0f, 0.0f, 0.0f);
  double3 vv = make_double3(0.0f, 0.0f, 0.0f);
  double ppx = 0;
  double ppy = 0;
  double ppz = 0;
  double vvx = 0;
  double vvy = 0;
  double vvz = 0;
  int i;

  // Calculate the total position and velocity of all particles
#pragma omp parallel for reduction(+:ppx,ppy,ppz,vvx,vvy,vvz)
  for (i = 0; i < n; ++i) {
    ppx += p[i].x * p[i].w;
    ppy += p[i].y * p[i].w;
    ppz += p[i].z * p[i].w;

    vvx += v[i].x * v[i].w;
    vvy += v[i].y * v[i].w;
    vvz += v[i].z * v[i].w;
  }

  // Calculate the center of mass position and velocity
  pp.x = ppx;
  pp.y = ppy;
  pp.z = ppz;

  vv.x = vvx;
  vv.y = vvy;
  vv.z = vvz;

  pp.x /= _M;
  pp.y /= _M;
  pp.z /= _M;
  vv.x /= _M;
  vv.y /= _M;
  vv.z /= _M;

  // Move particles to the center of mass and subtract the center of mass velocity
#pragma omp parallel for
  for (i = 0; i < n; ++i) {
    p[i].x -= pp.x;
    p[i].y -= pp.y;
    p[i].z -= pp.z;
    v[i].x -= vv.x;
    v[i].y -= vv.y;
    v[i].z -= vv.z;
  }
}

/**
 * Rescales the energy of the particles in the simulation.
 *
 * This function implements Algorithm 7.2 from Aarseth, 2003.
 * It rescales the potential energy and kinetic energy of the particles
 * to achieve a desired virial ratio and energy balance.
 *
 * @param n The number of particles.
 * @param p Pointer to the array of position vectors.
 * @param v Pointer to the array of velocity vectors.
 */
void RescaleEnergy(const int n, double4 *p, double4 *v) {
  // Aarseth, 2003, Algorithm 7.2.
  double Epot = Ep(n, p);
  double Ekin = Ek(n, v);
  double virialRatio = 0.5f;
  double Qv = sqrt(virialRatio * abs(Epot) / Ekin);
  scale4DArray(n, v, Qv);
  double beta = abs((1 - virialRatio) * Epot / (Epot + Ekin));

  scale4DArray(n, p, beta);
  scale4DArray(n, v, 1.0f / (sqrt(beta)));

  // After first scale Ekin is -0.5Epot but E0 != -0.25.
  // So just scale up or down as needed.
  Epot = Ep(n, p);
  beta = Epot / -0.5f;
  scale4DArray(n, p, beta);
  scale4DArray(n, v, 1.0f / sqrt(beta));
}

/**
 * Initializes the positions, velocities, and masses of the bodies.
 *
 * @param n The number of bodies.
 * @param p Pointer to the array of body positions.
 * @param v Pointer to the array of body velocities.
 * @param seed The seed for the random number generator (default is 0).
 */
void InitBodies(const int n, double4 *p, double4 *v, int seed = 0) {
  // Initialize masses equally
  InitPosVel(n, p, v, seed);

  // Translate bodies to move the center of mass on center of the coordinate
  // system
  Move2Center(n, p, v);

  // Rescale energy
  RescaleEnergy(n, p, v);
}

/**
 * @brief Computes the interactions between particles using the Barnes-Hut algorithm.
 *
 * This CUDA kernel calculates the forces between particles in a simulation using the Barnes-Hut algorithm.
 * It iterates over each particle and computes the forces between that particle and all other particles.
 * The forces are then used to update the velocities of the particles.
 *
 * @param n The number of particles in the simulation.
 * @param p Pointer to an array of double4 structures representing the positions and masses of the particles.
 * @param v Pointer to an array of double4 structures representing the velocities of the particles.
 * @param dt The time step for the simulation.
 */
__global__ void ComputeInteractions(const int n, double4 *p, double4 *v, const double dt) {
  double fx = 0.0f;
  double fy = 0.0f;
  double fz = 0.0f;
  auto i = blockDim.x * blockIdx.x + threadIdx.x;
  #pragma unroll
  for (int j = 0; j < n; ++j) {
      // compute distance pair
      auto dx = p[j].x - p[i].x;
      auto dy = p[j].y - p[i].y;
      auto dz = p[j].z - p[i].z;

      auto d = dx * dx + dy * dy + dz * dz + _SOFTENING*_SOFTENING;
      auto d_inv = 1.0f / sqrt(d);
      auto d_inv3 = d_inv * d_inv * d_inv;

      fx += d_inv3 * p[j].w * dx;
      fy += d_inv3 * p[j].w * dy;
      fz += d_inv3 * p[j].w * dz;  
  }
  v[i].x += fx * dt;
  v[i].y += fy * dt;
  v[i].z += fz * dt;
}

/**
 * @brief Updates the position of particles in a simulation.
 *
 * This CUDA kernel function updates the position of particles based on their velocity and the given time step.
 *
 * @param p - Pointer to an array of double4 structures representing the positions of particles.
 * @param v - Pointer to an array of double4 structures representing the velocities of particles.
 * @param dt - The time step used for the update.
 */
__global__ void UpdatePosition(double4 *p, double4 *v, const double dt) {
  auto i = blockDim.x * blockIdx.x + threadIdx.x;
  p[i].x += v[i].x * dt;
  p[i].y += v[i].y * dt;
  p[i].z += v[i].z * dt;
}


int main (int argc, char **argv) {
  int seed = 0;
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  if (argc == 3)
    seed = atoi(argv[2]);
  else  
    srand(seed);

  int n = atoi(argv[1]);
  const double dt = 0.01f; 
  
  double4 *h_p, *h_v;

  // Set device specified as argument to the program
  setDeviceByName(argv[3]);

  // Allocate pinned memory
  cudaMallocHost(&h_p, n * sizeof(double4));
  cudaMallocHost(&h_v, n * sizeof(double4));

  // Init Bodies
  TIMERSTART(init)
  InitBodies(n, h_p, h_v, seed);
#ifdef MONITOR_ENERGY
  {
  double ek = Ek(n, h_v);
  double ep = Ep(n, h_p);
  std::cout << "Etot: " <<ek+ep <<std::endl;
  }
#endif
  TIMERSTOP(init)

  // Allocate memory on the device
  double4 *d_p, *d_v;
  cudaMalloc(&d_p, n * sizeof(double4));
  cudaMalloc(&d_v, n * sizeof(double4));

  dim3 blocks, threadsPerBlock;
  if (n < MAX_THREADS_PER_BLOCK) {
    blocks = dim3(1);
    threadsPerBlock = dim3(n);
  } else {
    // assuming that n is a power of two
    blocks = dim3(n / MAX_THREADS_PER_BLOCK);
    threadsPerBlock = dim3(MAX_THREADS_PER_BLOCK);
  }

  TIMERSTART(total)
  cudaMemcpy(d_p, h_p, n * sizeof(double4), cudaMemcpyHostToDevice);
  cudaMemcpy(d_v, h_v, n * sizeof(double4), cudaMemcpyHostToDevice);
  
  TIMERSTART(simulation)
  for (double t = 0; t < 100; t+= dt) {
    ComputeInteractions<<<blocks, threadsPerBlock>>>(n, d_p, d_v, dt );
    UpdatePosition<<<blocks, threadsPerBlock>>>(d_p, d_v, dt);
    cudaMemcpy(h_p, d_p, n * sizeof(double4), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_v, d_v, n * sizeof(double4), cudaMemcpyDeviceToHost);
#ifdef MONITOR_ENERGY
    {
      double ek = Ek(n, h_v);
      double ep = Ep(n, h_p);
      std::cout << "Etot: " <<ek+ep <<std::endl;
    }
#endif
  }
  cudaDeviceSynchronize();
  TIMERSTOP(simulation)

  //cudaMemcpy(h_p, d_p, n * sizeof(double4), cudaMemcpyDeviceToHost);
  //cudaMemcpy(h_v, d_v, n * sizeof(double4), cudaMemcpyDeviceToHost);
  TIMERSTOP(total)

#ifdef MONITOR_ENERGY
  {
    double ek = Ek(n, h_v);
    double ep = Ep(n, h_p);
    std::cout << "Etot: " <<ek+ep <<std::endl;
  }
#endif
  
  // Free pinned memory
  cudaFreeHost(h_p);
  cudaFreeHost(h_v);

  // Free memory on device
  cudaFree(d_p);
  cudaFree(d_v);

  return 0;
}
