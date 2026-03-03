#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

#include "gpu_manhattan_cuda.h"

#ifndef CELL_CAP
#define CELL_CAP 16
#endif

static inline void ck(cudaError_t e, const char* msg) {
    if (e != cudaSuccess) {
        fprintf(stderr, "CUDA error (%s): %s\n", msg, cudaGetErrorString(e));
        std::exit(1);
    }
}

static int g_inited = 0;

__device__ __forceinline__ int wrap(int a, int nc) {
    if (a < 0) return a + nc;
    if (a >= nc) return a - nc;
    return a;
}
__device__ __forceinline__ int cell_index(int x, int y, int z, int nc) {
    return x + nc * (y + nc * z);
}

/**
 * Owner-computes-i (no atomics for forces): thread i accumulates force on i.
 * Potential energy: we sum full pair energy into Ui (double counts pairs),
 * then multiply by 0.5 on the host to match standard U.
 */
__global__ void manhattan_forces_U_owner_i(
    const float* __restrict__ x,
    const float* __restrict__ y,
    const float* __restrict__ z,
    float* __restrict__ fx,
    float* __restrict__ fy,
    float* __restrict__ fz,
    float* __restrict__ Ui,
    const int* __restrict__ counts,
    const int* __restrict__ cells,
    unsigned N,
    int nc,
    float L,
    float rc2,
    float sigma,
    float epsilon,
    float cell_size
){
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    float xi = x[i], yi = y[i], zi = z[i];

    int cx = (int)floorf(xi / cell_size);
    int cy = (int)floorf(yi / cell_size);
    int cz = (int)floorf(zi / cell_size);
    cx = wrap(cx, nc); cy = wrap(cy, nc); cz = wrap(cz, nc);

    int self = cell_index(cx, cy, cz, nc);

    int neigh[7];
    neigh[0] = self;
    neigh[1] = cell_index(wrap(cx+1,nc), cy, cz, nc);
    neigh[2] = cell_index(wrap(cx-1,nc), cy, cz, nc);
    neigh[3] = cell_index(cx, wrap(cy+1,nc), cz, nc);
    neigh[4] = cell_index(cx, wrap(cy-1,nc), cz, nc);
    neigh[5] = cell_index(cx, cy, wrap(cz+1,nc), nc);
    neigh[6] = cell_index(cx, cy, wrap(cz-1,nc), nc);

    float fxi = 0.f, fyi = 0.f, fzi = 0.f;
    float ui  = 0.f;

    for (int n=0; n<7; n++) {
        int c = neigh[n];
        int cnt = counts[c];
        int base = c * CELL_CAP;

        for (int k=0; k<cnt; k++) {
            int j = cells[base + k];
            if (j < 0 || (unsigned)j == i) continue;

            float dx = xi - x[j];
            float dy = yi - y[j];
            float dz = zi - z[j];

            // minimum image
            dx -= L * nearbyintf(dx / L);
            dy -= L * nearbyintf(dy / L);
            dz -= L * nearbyintf(dz / L);

            float r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > rc2 || r2 < 1e-20f) continue;

            float inv_r2  = 1.0f / r2;
            float sig2_r2 = (sigma*sigma) * inv_r2;
            float sr6     = sig2_r2 * sig2_r2 * sig2_r2;
            float sr12    = sr6 * sr6;

            // Potential (pair)
            ui += 4.0f * epsilon * (sr12 - sr6);

            // Force: F = 24ε(2sr12 - sr6) * (1/r^2) * r_vec
            float F_over_r = 24.0f * epsilon * (2.0f*sr12 - sr6) * inv_r2;

            fxi += F_over_r * dx;
            fyi += F_over_r * dy;
            fzi += F_over_r * dz;
        }
    }

    fx[i] = fxi;
    fy[i] = fyi;
    fz[i] = fzi;
    Ui[i] = ui;
}

int gpu_manhattan_cuda_init(void) {
    if (g_inited) return 1;
    int ndev = 0;
    cudaError_t e = cudaGetDeviceCount(&ndev);
    if (e != cudaSuccess || ndev <= 0) return 0;
    g_inited = 1;
    return 1;
}

void gpu_manhattan_cuda_cleanup(void) {
    g_inited = 0;
}

double gpu_manhattan_cuda_forces_and_U(Particle *p, const SimParams *sp, const CellList *cl) {
    if (!g_inited) return 0.0;

    const unsigned N = (unsigned)sp->N;

    // host SoA positions
    std::vector<float> hx(N), hy(N), hz(N);
    for (unsigned i=0;i<N;i++){
        hx[i] = (float)p[i].x;
        hy[i] = (float)p[i].y;
        hz[i] = (float)p[i].z;
    }

    float *dx=nullptr, *dy=nullptr, *dz=nullptr;
    float *dfx=nullptr, *dfy=nullptr, *dfz=nullptr;
    float *dUi=nullptr;
    int *dcounts=nullptr, *dcells=nullptr;

    const int ncell = cl->ncell;
    const size_t cells_bytes = (size_t)ncell * CELL_CAP * sizeof(int);

    ck(cudaMalloc(&dx,  sizeof(float)*N), "malloc x");
    ck(cudaMalloc(&dy,  sizeof(float)*N), "malloc y");
    ck(cudaMalloc(&dz,  sizeof(float)*N), "malloc z");
    ck(cudaMalloc(&dfx, sizeof(float)*N), "malloc fx");
    ck(cudaMalloc(&dfy, sizeof(float)*N), "malloc fy");
    ck(cudaMalloc(&dfz, sizeof(float)*N), "malloc fz");
    ck(cudaMalloc(&dUi, sizeof(float)*N), "malloc Ui");

    ck(cudaMalloc(&dcounts, sizeof(int)*ncell), "malloc counts");
    ck(cudaMalloc(&dcells,  cells_bytes),       "malloc cells");

    ck(cudaMemcpy(dx, hx.data(), sizeof(float)*N, cudaMemcpyHostToDevice), "cpy x");
    ck(cudaMemcpy(dy, hy.data(), sizeof(float)*N, cudaMemcpyHostToDevice), "cpy y");
    ck(cudaMemcpy(dz, hz.data(), sizeof(float)*N, cudaMemcpyHostToDevice), "cpy z");
    ck(cudaMemcpy(dcounts, cl->counts, sizeof(int)*ncell, cudaMemcpyHostToDevice), "cpy counts");
    ck(cudaMemcpy(dcells,  cl->cells,  cells_bytes, cudaMemcpyHostToDevice),       "cpy cells");

    const int threads = 256;
    const int blocks  = (int)((N + threads - 1) / threads);

    manhattan_forces_U_owner_i<<<blocks, threads>>>(
        dx, dy, dz, dfx, dfy, dfz, dUi,
        dcounts, dcells,
        N,
        cl->nc,
        (float)sp->L,
        (float)sp->rc2,
        (float)sp->sigma,
        (float)sp->epsilon,
        (float)cl->cell_size
    );
    ck(cudaGetLastError(), "kernel launch");
    ck(cudaDeviceSynchronize(), "sync");

    std::vector<float> hfx(N), hfy(N), hfz(N), hUi(N);
    ck(cudaMemcpy(hfx.data(), dfx, sizeof(float)*N, cudaMemcpyDeviceToHost), "cpy fx back");
    ck(cudaMemcpy(hfy.data(), dfy, sizeof(float)*N, cudaMemcpyDeviceToHost), "cpy fy back");
    ck(cudaMemcpy(hfz.data(), dfz, sizeof(float)*N, cudaMemcpyDeviceToHost), "cpy fz back");
    ck(cudaMemcpy(hUi.data(), dUi, sizeof(float)*N, cudaMemcpyDeviceToHost), "cpy Ui back");

    for (unsigned i=0;i<N;i++){
        p[i].fx = (double)hfx[i];
        p[i].fy = (double)hfy[i];
        p[i].fz = (double)hfz[i];
    }

    // Ui contains double-counted pair energy (i sums over all j).
    // Multiply by 0.5 to match conventional total U.
    double U = 0.0;
    for (unsigned i=0;i<N;i++) U += (double)hUi[i];
    U *= 0.5;

    cudaFree(dx); cudaFree(dy); cudaFree(dz);
    cudaFree(dfx); cudaFree(dfy); cudaFree(dfz);
    cudaFree(dUi);
    cudaFree(dcounts); cudaFree(dcells);

    return U;
}