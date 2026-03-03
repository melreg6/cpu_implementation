#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <cuda_runtime.h>
#include <cub/cub.cuh>

#include "gpu_manhattan.h"

#ifndef CUDA_CHECK
#define CUDA_CHECK(call) do {                                \
  cudaError_t err = (call);                                  \
  if (err != cudaSuccess) {                                  \
    fprintf(stderr, "CUDA error %s:%d: %s\n",                \
            __FILE__, __LINE__, cudaGetErrorString(err));    \
    std::exit(1);                                            \
  }                                                          \
} while (0)
#endif

// -------------------------
// helpers (device)
// -------------------------
__device__ __forceinline__ double wrap_pbc(double x, double L) {
  // map to [0, L)
  x = fmod(x, L);
  if (x < 0) x += L;
  return x;
}

__device__ __forceinline__ double min_image(double dx, double L) {
  // minimum image convention
  if (dx >  0.5 * L) dx -= L;
  if (dx < -0.5 * L) dx += L;
  return dx;
}

__device__ __forceinline__ int3 cell_coord(double x, double y, double z, double L, int nc) {
  x = wrap_pbc(x, L);
  y = wrap_pbc(y, L);
  z = wrap_pbc(z, L);
  double cs = L / (double)nc;
  int cx = (int)floor(x / cs);
  int cy = (int)floor(y / cs);
  int cz = (int)floor(z / cs);
  if (cx >= nc) cx = nc - 1;
  if (cy >= nc) cy = nc - 1;
  if (cz >= nc) cz = nc - 1;
  return make_int3(cx, cy, cz);
}

__device__ __forceinline__ int cell_id_from_coord(int cx, int cy, int cz, int nc) {
  cx = (cx % nc + nc) % nc;
  cy = (cy % nc + nc) % nc;
  cz = (cz % nc + nc) % nc;
  return (cz * nc + cy) * nc + cx;
}

// -------------------------
// kernels
// -------------------------
__global__ void k_zero_forces(Particle* p, size_t N) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    p[i].fx = 0.0;
    p[i].fy = 0.0;
    p[i].fz = 0.0;
  }
}

__global__ void k_compute_cell_ids(const Particle* p, int* cellId, int* pid,
                                   size_t N, double L, int nc) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    int3 c = cell_coord(p[i].x, p[i].y, p[i].z, L, nc);
    int id = cell_id_from_coord(c.x, c.y, c.z, nc);
    cellId[i] = id;
    pid[i] = (int)i;
  }
}

__global__ void k_init_cell_ranges(int* cellStart, int* cellEnd, int ncell) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < ncell) {
    cellStart[i] = -1;
    cellEnd[i]   = -1;
  }
}

__global__ void k_build_cell_ranges(const int* sortedCellId, int* cellStart, int* cellEnd,
                                    size_t N) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= N) return;

  int cid = sortedCellId[i];

  if (i == 0 || sortedCellId[i - 1] != cid) {
    cellStart[cid] = (int)i;
  }
  if (i == N - 1 || sortedCellId[i + 1] != cid) {
    cellEnd[cid] = (int)i + 1;
  }
}

__global__ void k_forces_manhattan(
    Particle* p,
    const int* sortedPid,
    const int* cellStart,
    const int* cellEnd,
    size_t N,
    double L,
    int nc,
    double rc2,
    double sigma,
    double epsilon,
    double* U_partial
) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= N) return;

  double xi = p[i].x, yi = p[i].y, zi = p[i].z;
  int3 ci = cell_coord(xi, yi, zi, L, nc);

  // Manhattan radius-1: self + 6 face neighbors
  const int dxs[7] = {0, 1, -1, 0, 0, 0, 0};
  const int dys[7] = {0, 0,  0, 1,-1, 0, 0};
  const int dzs[7] = {0, 0,  0, 0, 0, 1,-1};

  double fx = 0.0, fy = 0.0, fz = 0.0;
  double Ui = 0.0;

  #pragma unroll
  for (int nb = 0; nb < 7; nb++) {
    int cid = cell_id_from_coord(ci.x + dxs[nb], ci.y + dys[nb], ci.z + dzs[nb], nc);

    int start = cellStart[cid];
    int end   = cellEnd[cid];
    if (start < 0 || end < 0) continue;

    for (int s = start; s < end; s++) {
      int j = sortedPid[s];
      if (j == (int)i) continue;

      double dx = min_image(p[j].x - xi, L);
      double dy = min_image(p[j].y - yi, L);
      double dz = min_image(p[j].z - zi, L);

      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 <= rc2 && r2 > 1e-20) {
        double inv_r2 = 1.0 / r2;
        double sig2 = sigma * sigma;
        double sr2 = sig2 * inv_r2;
        double sr6 = sr2 * sr2 * sr2;
        double sr12 = sr6 * sr6;

        Ui += 4.0 * epsilon * (sr12 - sr6);

        double coeff = 24.0 * epsilon * (2.0*sr12 - sr6) * inv_r2;
        fx += coeff * dx;
        fy += coeff * dy;
        fz += coeff * dz;
      }
    }
  }

  p[i].fx = fx;
  p[i].fy = fy;
  p[i].fz = fz;
  U_partial[i] = Ui;
}

// -------------------------
// host entry point
// -------------------------
double md_compute_forces_cell_manhattan_gpu(Particle* p_h, const SimParams* sp_h) {
  const size_t N = sp_h->N;
  const double L = sp_h->L;

  int nc = (int)floor(L / sp_h->rc);
  if (nc < 1) nc = 1;
  int ncell = nc * nc * nc;

  Particle* p_d = nullptr;
  int *cellId_d = nullptr, *pid_d = nullptr;
  int *cellId_out_d = nullptr, *pid_out_d = nullptr;
  int *cellStart_d = nullptr, *cellEnd_d = nullptr;
  double* U_partial_d = nullptr;

  CUDA_CHECK(cudaMalloc(&p_d, N * sizeof(Particle)));
  CUDA_CHECK(cudaMemcpy(p_d, p_h, N * sizeof(Particle), cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMalloc(&cellId_d, N * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&pid_d,    N * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&cellStart_d, ncell * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&cellEnd_d,   ncell * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&U_partial_d, N * sizeof(double)));

  const int block = 256;
  int gridN = (int)((N + block - 1) / block);
  int gridCell = (ncell + block - 1) / block;

  k_zero_forces<<<gridN, block>>>(p_d, N);
  CUDA_CHECK(cudaGetLastError());

  k_compute_cell_ids<<<gridN, block>>>(p_d, cellId_d, pid_d, N, L, nc);
  CUDA_CHECK(cudaGetLastError());

  // ---- CUB radix sort pairs (cellId, pid)
  CUDA_CHECK(cudaMalloc(&cellId_out_d, N * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&pid_out_d,    N * sizeof(int)));

  void* d_temp = nullptr;
  size_t temp_bytes = 0;

  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(
      d_temp, temp_bytes,
      cellId_d, cellId_out_d,
      pid_d, pid_out_d,
      (int)N
  ));

  CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));

  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(
      d_temp, temp_bytes,
      cellId_d, cellId_out_d,
      pid_d, pid_out_d,
      (int)N
  ));

  CUDA_CHECK(cudaFree(d_temp));
  CUDA_CHECK(cudaFree(cellId_d));
  CUDA_CHECK(cudaFree(pid_d));
  cellId_d = cellId_out_d;
  pid_d = pid_out_d;

  k_init_cell_ranges<<<gridCell, block>>>(cellStart_d, cellEnd_d, ncell);
  CUDA_CHECK(cudaGetLastError());

  k_build_cell_ranges<<<gridN, block>>>(cellId_d, cellStart_d, cellEnd_d, N);
  CUDA_CHECK(cudaGetLastError());

  k_forces_manhattan<<<gridN, block>>>(
      p_d, pid_d, cellStart_d, cellEnd_d, N,
      L, nc, sp_h->rc2, sp_h->sigma, sp_h->epsilon, U_partial_d
  );
  CUDA_CHECK(cudaGetLastError());

  std::vector<double> U_partial_h(N);
  CUDA_CHECK(cudaMemcpy(U_partial_h.data(), U_partial_d, N*sizeof(double), cudaMemcpyDeviceToHost));

  double U = 0.0;
  for (size_t i = 0; i < N; i++) U += U_partial_h[i];
  U *= 0.5; // double-count correction

  CUDA_CHECK(cudaMemcpy(p_h, p_d, N * sizeof(Particle), cudaMemcpyDeviceToHost));

  CUDA_CHECK(cudaFree(p_d));
  CUDA_CHECK(cudaFree(cellId_d));
  CUDA_CHECK(cudaFree(pid_d));
  CUDA_CHECK(cudaFree(cellStart_d));
  CUDA_CHECK(cudaFree(cellEnd_d));
  CUDA_CHECK(cudaFree(U_partial_d));

  return U;
}