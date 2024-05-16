#include <hip/hip_runtime.h>
#include "utilities/hip_utils.h"
#include "utilities/time_utils.h"
#include "utilities/nbody_helpers.h"

#define MAX_THREADS_PER_BLOCK 32

__global__ void ComputeInteractions(const int n, float *m, float *p, float *v, const float dt) {
  float fx = 0.0f;
  float fy = 0.0f;
  float fz = 0.0f;
  auto i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x * 3;
  for (int j = 0; j < 3*n; j+=3) {
      auto m2_id = j / 3;
      // compute distance pair
      auto dx = p[j] - p[i];
      auto dy = p[j + 1] - p[i + 1];
      auto dz = p[j + 2] - p[i + 2];

      auto d = dx * dx + dy * dy + dz * dz + _SOFTENING*_SOFTENING;
      auto d_inv = 1.0f / sqrtf(d);
      auto d_inv3 = d_inv * d_inv * d_inv;

      fx += d_inv3 * m[m2_id] * dx;
      fy += d_inv3 * m[m2_id] * dy;
      fz += d_inv3 * m[m2_id] * dz;  
  }
  v[i] += fx * dt;
  v[i + 1] += fy * dt;
  v[i + 2] += fz * dt;
}

__global__ void UpdatePosition(float *p, float *v, const float dt) {
  auto i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x * 3;
  p[i] += v[i] * dt;
  p[i+1] += v[i+1] * dt;
  p[i+2] += v[i+2] * dt;
}

int main (int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  srand(0);
  
  int n = atoi(argv[1]);
  const float dt = 0.01f; 

  float *h_m, *h_p, *h_v;

  // Allocate pinned memory
  HIP_CHECK(hipHostMalloc(&h_m, n * sizeof(float)))
  HIP_CHECK(hipHostMalloc(&h_p, 3*n * sizeof(float)))
  HIP_CHECK(hipHostMalloc(&h_v, 3*n * sizeof(float)))

  // Init Bodies
  InitAos(n, h_m, h_p, h_v);

  // Allocate memory on the device
  float *d_m, *d_p, *d_v;
  HIP_CHECK(hipMalloc(&d_m, n * sizeof(float)))
  HIP_CHECK(hipMalloc(&d_p, 3*n * sizeof(float)))
  HIP_CHECK(hipMalloc(&d_v, 3*n * sizeof(float)))
  
  HIP_CHECK(hipMemcpy(d_m, h_m, n * sizeof(float), hipMemcpyHostToDevice))
  HIP_CHECK(hipMemcpy(d_p, h_p, 3*n * sizeof(float), hipMemcpyHostToDevice))
  HIP_CHECK(hipMemcpy(d_v, h_v, 3*n * sizeof(float), hipMemcpyHostToDevice))
  
  dim3 blocks, threadsPerBlock;

  if (n < MAX_THREADS_PER_BLOCK) {
    blocks = dim3(1);
    threadsPerBlock = dim3(n);
  } else {
    // assuming that n is a power of two
    blocks = dim3(n / MAX_THREADS_PER_BLOCK);
    threadsPerBlock = dim3(MAX_THREADS_PER_BLOCK);
  }
   
  TIMERSTART(simulation)
  for (float t = 0; t < 0.1; t+= dt) {
    hipLaunchKernelGGL(ComputeInteractions /* compute kernel*/, 
                        blocks, threadsPerBlock, 0/*dynamic shared*/, 0/*stream*/, 
                        n, d_m, d_p, d_v, dt /* arguments to the compute kernel */);
    HIP_CHECK(hipDeviceSynchronize())
    hipLaunchKernelGGL(UpdatePosition,
                        blocks, threadsPerBlock, 0/*dynamic shared*/, 0/*stream*/, 
                        d_p, d_v, dt /* arguments to the compute kernel */);
    HIP_CHECK(hipDeviceSynchronize())
  }
  TIMERSTOP(simulation)

  HIP_CHECK(hipMemcpy(h_m, d_m, n * sizeof(float), hipMemcpyDeviceToHost))
  HIP_CHECK(hipMemcpy(h_p, d_p, 3*n * sizeof(float), hipMemcpyDeviceToHost))
  HIP_CHECK(hipMemcpy(h_v, d_v, 3*n * sizeof(float), hipMemcpyDeviceToHost))

  float ek = Ek(n, h_m, h_v);
  float ep = Ep(n, h_m, h_p);
  std::cout << "Etot: " <<ek+ep <<std::endl;
  
  // Free pinned memory
  HIP_CHECK(hipHostFree(h_m))
  HIP_CHECK(hipHostFree(h_p))
  HIP_CHECK(hipHostFree(h_v))

  // Free memory on device
  HIP_CHECK(hipFree(d_m))
  HIP_CHECK(hipFree(d_p))
  HIP_CHECK(hipFree(d_v))

  return 0;
}