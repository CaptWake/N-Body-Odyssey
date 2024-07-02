#include "utilities/time_utils.h"
#include "utilities/nbody_helpers.h"
#include <iostream>
#include <algorithm>
#include <omp.h>

#define MAX_THREADS_PER_BLOCK 1024

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

float inline xoxyi_rand(unsigned int *seed){
  return (float)rand_r(seed) / (float) RAND_MAX;
}

float inline Ep(const int n, float4 *p) {
  float Epot = 0.0f;
#pragma omp parallel for reduction(+:Epot)
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      float D, x, y, z;
      x = p[i].x - p[j].x;
      y = p[i].y - p[j].y;
      z = p[i].z - p[j].z;
      D = sqrtf(x * x + y * y + z * z);
      Epot += -1.0f * p[i].w * p[j].w / D;
    }
  }
  return Epot;
}

float inline Ek(const int n, float4 *v) {
  float Ekin = 0.0;
#pragma omp parallel for reduction(+:Ekin)
  for (int i = 0; i < n; ++i) {
    Ekin += 0.5f * v[i].w * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z);
  }
  return Ekin;
}

static inline void scale4DArray(const int n, float4 *m, const float scale) {
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    m[i].x *= scale;
    m[i].y *= scale;
    m[i].z *= scale;
  }
}

void InitPosVel(const int n, float4 *p, float4 *v, int seed) {
#pragma omp parallel
 { 
  unsigned int seedT = omp_get_thread_num() + seed * omp_get_num_threads();	 
  #pragma omp for
    for (int i = 0; i < n; ++i) {
      float R, X, Y;
      float mi = _M / n;
      R = xoxyi_rand(&seedT);
      X = acosf(1.0f - 2.0f * xoxyi_rand(&seedT));
      Y = xoxyi_rand(&seedT) * 2.0f * _PI;
      // https://www.researchgate.net/figure/Figure-A1-Spherical-coordinates_fig8_284609648
      p[i] = make_float4(R * sinf(X) * cosf(Y), R * sinf(X) * sinf(Y), R * cosf(X), mi);
      v[i] = make_float4(1.0f - 2.0f * xoxyi_rand(&seedT), 1.0f - 2.0f * xoxyi_rand(&seedT), 1.0f - 2.0f * xoxyi_rand(&seedT), mi);
    }
  }
}

void Move2Center(const int n, float4 *p, float4 *v) {
  float3 pp = make_float3(0.0f, 0.0f, 0.0f);
  float3 vv = make_float3(0.0f, 0.0f, 0.0f);
  float ppx = 0;
  float ppy = 0;
  float ppz = 0;
  float vvx = 0;
  float vvy = 0;
  float vvz = 0;
  int i;
#pragma omp parallel for reduction(+:ppx,ppy,ppz,vvx,vvy,vvz)
  for (i = 0; i < n; ++i) {
    ppx += p[i].x * p[i].w;
    ppy += p[i].y * p[i].w;
    ppz += p[i].z * p[i].w;

    vvx += v[i].x * v[i].w;
    vvy += v[i].y * v[i].w;
    vvz += v[i].z * v[i].w;
  }

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

void RescaleEnergy(const int n, float4 *p, float4 *v) {
  // Aarseth, 2003, Algorithm 7.2.
  float Epot = Ep(n, p);
  float Ekin = Ek(n, v);
  float virialRatio = 0.5f;
  float Qv = sqrtf(virialRatio * fabsf(Epot) / Ekin);
  scale4DArray(n, v, Qv);
  float beta = fabsf((1 - virialRatio) * Epot / (Epot + Ekin));

  scale4DArray(n, p, beta);
  scale4DArray(n, v, 1.0f / (sqrtf(beta)));

  // After first scale Ekin is -0.5Epot but E0 != -0.25.
  // So just scale up or down as needed.
  Epot = Ep(n, p);
  beta = Epot / -0.5f;
  scale4DArray(n, p, beta);
  scale4DArray(n, v, 1.0f / sqrtf(beta));
}

void InitBodies(const int n, float4 *p, float4 *v, int seed = 0) {
  // Initialize masses equally
  InitPosVel(n, p, v, seed);

  // Translate bodies to move the center of mass on center of the coordinate
  // system
  Move2Center(n, p, v);

  // Rescale energy
  RescaleEnergy(n, p, v);
}

__global__ void ComputeInteractions(const int n, float4 *p, float4 *v, const float dt) {
  float fx = 0.0f;
  float fy = 0.0f;
  float fz = 0.0f;
  auto i = blockDim.x * blockIdx.x + threadIdx.x;
  #pragma unroll
  for (int j = 0; j < n; ++j) {
      // compute distance pair
      auto dx = p[j].x - p[i].x;
      auto dy = p[j].y - p[i].y;
      auto dz = p[j].z - p[i].z;

      auto d = dx * dx + dy * dy + dz * dz + _SOFTENING*_SOFTENING;
      auto d_inv = 1.0f / sqrtf(d);
      auto d_inv3 = d_inv * d_inv * d_inv;

      fx += d_inv3 * p[j].w * dx;
      fy += d_inv3 * p[j].w * dy;
      fz += d_inv3 * p[j].w * dz;  
  }
  v[i].x += fx * dt;
  v[i].y += fy * dt;
  v[i].z += fz * dt;
}

__global__ void UpdatePosition(float4 *p, float4 *v, const float dt) {
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
  const float dt = 0.01f; 
  
  float4 *h_p, *h_v;

  setDeviceByName(argv[3]);

  // Allocate pinned memory
  cudaMallocHost(&h_p, n * sizeof(float4));
  cudaMallocHost(&h_v, n * sizeof(float4));

  // Init Bodies
  TIMERSTART(init)
  InitBodies(n, h_p, h_v, seed);
#ifdef MONITOR_ENERGY
  {
  float ek = Ek(n, h_v);
  float ep = Ep(n, h_p);
  std::cout << "Etot: " <<ek+ep <<std::endl;
  }
#endif
  TIMERSTOP(init)

  // Allocate memory on the device
  float4 *d_p, *d_v;
  cudaMalloc(&d_p, n * sizeof(float4));
  cudaMalloc(&d_v, n * sizeof(float4));

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
  cudaMemcpy(d_p, h_p, n * sizeof(float4), cudaMemcpyHostToDevice);
  cudaMemcpy(d_v, h_v, n * sizeof(float4), cudaMemcpyHostToDevice);
  
  TIMERSTART(simulation)
  for (float t = 0; t < 100; t+= dt) {
    ComputeInteractions<<<blocks, threadsPerBlock>>>(n, d_p, d_v, dt );
    UpdatePosition<<<blocks, threadsPerBlock>>>(d_p, d_v, dt);
    cudaMemcpy(h_p, d_p, n * sizeof(float4), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_v, d_v, n * sizeof(float4), cudaMemcpyDeviceToHost);
#ifdef MONITOR_ENERGY
    {
      float ek = Ek(n, h_v);
      float ep = Ep(n, h_p);
      std::cout << "Etot: " <<ek+ep <<std::endl;
    }
#endif
  }
  cudaDeviceSynchronize();
  TIMERSTOP(simulation)

  //cudaMemcpy(h_p, d_p, n * sizeof(float4), cudaMemcpyDeviceToHost);
  //cudaMemcpy(h_v, d_v, n * sizeof(float4), cudaMemcpyDeviceToHost);
  TIMERSTOP(total)
#ifdef MONITOR_ENERGY
  {
  float ek = Ek(n, h_v);
  float ep = Ep(n, h_p);
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
