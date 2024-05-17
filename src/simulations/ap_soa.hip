#include <hip/hip_runtime.h>
#include "utilities/hip_utils.h"
#include "utilities/time_utils.h"
#include "utilities/nbody_helpers.h"

#define MAX_THREADS_PER_BLOCK 32

float inline Ep(const int n, float4 *p) {
  float Epot = 0.0f;
  float D, x, y, z;
  uint64_t i, j;
  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
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
  int i;
  for (i = 0; i < n; ++i) {
    Ekin += 0.5f * v[i].w * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z);
  }
  return Ekin;
}

static inline void scale4DArray(const int n, float4 *m, const float scale) {
  for (int i = 0; i < n; ++i) {
    m[i].x *= scale;
    m[i].y *= scale;
    m[i].z *= scale;
  }
}

void InitPosVel(const int n, float4 *p, float4 *v) {
  float R, X, Y;
  float mi = _M / n;
  for (int i = 0; i < n; ++i) {
    R = frand();
    X = acosf(1.0f - 2.0f * frand());
    Y = frand() * 2.0f * _PI;
    // https://www.researchgate.net/figure/Figure-A1-Spherical-coordinates_fig8_284609648
    p[i] = make_float4(R * sinf(X) * cosf(Y), R * sinf(X) * sinf(Y), R * cosf(X), mi);
    v[i] = make_float4(1.0f - 2.0f * frand(), 1.0f - 2.0f * frand(), 1.0f - 2.0f * frand(), mi);
  }
}

void Move2Center(const int n, float4 *p, float4 *v) {
  float3 pp = make_float3(0.0f, 0.0f, 0.0f);
  float3 vv = make_float3(0.0f, 0.0f, 0.0f);
  int i;
  for (i = 0; i < n; ++i) {
    pp.x += p[i].x * p[i].w;
    pp.y += p[i].y * p[i].w;
    pp.z += p[i].z * p[i].w;

    vv.x += v[i].x * v[i].w;
    vv.y += v[i].y * v[i].w;
    vv.z += v[i].z * v[i].w;
  }

  pp.x /= _M;
  pp.y /= _M;
  pp.z /= _M;
  vv.x /= _M;
  vv.y /= _M;
  vv.z /= _M;

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

void InitBodies(const int n, float4 *p, float4 *v) {
  // Initialize masses equally
  InitPosVel(n, p, v);

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
  auto i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
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
  auto i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  p[i].x += v[i].x * dt;
  p[i].y += v[i].y * dt;
  p[i].z += v[i].z * dt;
}


int main (int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Must specify the number of bodies" << std::endl;
    exit(1);
  }
  srand(0);
  
  int n = atoi(argv[1]);
  const float dt = 0.01f; 
  
  float4 *h_p, *h_v;

  // Allocate pinned memory
  HIP_CHECK(hipHostMalloc(&h_p, n * sizeof(float4)))
  HIP_CHECK(hipHostMalloc(&h_v, n * sizeof(float4)))

  // Init Bodies
  InitBodies(n, h_p, h_v);

  // Allocate memory on the device
  float4 *d_p, *d_v;
  HIP_CHECK(hipMalloc(&d_p, n * sizeof(float4)))
  HIP_CHECK(hipMalloc(&d_v, n * sizeof(float4)))
  
  HIP_CHECK(hipMemcpy(d_p, h_p, n * sizeof(float4), hipMemcpyHostToDevice))
  HIP_CHECK(hipMemcpy(d_v, h_v, n * sizeof(float4), hipMemcpyHostToDevice))
  
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
                        n, d_p, d_v, dt /* arguments to the compute kernel */);
    HIP_CHECK(hipDeviceSynchronize())
    hipLaunchKernelGGL(UpdatePosition,
                        blocks, threadsPerBlock, 0/*dynamic shared*/, 0/*stream*/, 
                        d_p, d_v, dt /* arguments to the compute kernel */);
    HIP_CHECK(hipDeviceSynchronize())
  }
  TIMERSTOP(simulation)

  HIP_CHECK(hipMemcpy(h_p, d_p, n * sizeof(float4), hipMemcpyDeviceToHost))
  HIP_CHECK(hipMemcpy(h_v, d_v, n * sizeof(float4), hipMemcpyDeviceToHost))

  float ek = Ek(n, h_v);
  float ep = Ep(n, h_p);
  std::cout << "Etot: " <<ek+ep <<std::endl;
  
  // Free pinned memory
  HIP_CHECK(hipHostFree(h_p))
  HIP_CHECK(hipHostFree(h_v))

  // Free memory on device
  HIP_CHECK(hipFree(d_p))
  HIP_CHECK(hipFree(d_v))

  return 0;
}