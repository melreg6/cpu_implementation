#include "cell_manhattan_pthread.h"

#include <pthread.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef CELL_CAP
// MUST match your CellList fixed capacity per cell.
// Your existing code uses 16 in multiple places.
#define CELL_CAP 16
#endif

static inline int cell_index(int ix, int iy, int iz, int nc) {
    return ix + nc * (iy + nc * iz);
}

static inline int wrap(int a, int nc) {
    if (a < 0) return a + nc;
    if (a >= nc) return a - nc;
    return a;
}

// 7-cell Manhattan neighborhood (radius 1): self + ±x ±y ±z
static inline void manhattan_neighbors_7(const CellList *cl, int cell, int out7[7]) {
    const int nc = cl->nc;

    int z = cell / (nc * nc);
    int rem = cell % (nc * nc);
    int y = rem / nc;
    int x = rem % nc;

    out7[0] = cell;
    out7[1] = cell_index(wrap(x + 1, nc), y, z, nc);
    out7[2] = cell_index(wrap(x - 1, nc), y, z, nc);
    out7[3] = cell_index(x, wrap(y + 1, nc), z, nc);
    out7[4] = cell_index(x, wrap(y - 1, nc), z, nc);
    out7[5] = cell_index(x, y, wrap(z + 1, nc), nc);
    out7[6] = cell_index(x, y, wrap(z - 1, nc), nc);
}

typedef struct {
    int cell_start;
    int cell_end;

    Particle *p;
    const SimParams *sp;
    const CellList *cl;

    double *local_fx;
    double *local_fy;
    double *local_fz;

    double local_U;
} manh_thread_data_t;

static void* thread_compute_force_manhattan(void *arg) {
    manh_thread_data_t *ctx = (manh_thread_data_t*)arg;

    Particle *p = ctx->p;
    const SimParams *sp = ctx->sp;
    const CellList *cl = ctx->cl;

    const double L   = sp->L;
    const double rc2 = sp->rc2;
    const double sig = sp->sigma;
    const double eps = sp->epsilon;

    const double tiny2 = 1e-30;
    double U = 0.0;

    int neigh[7];

    for (int cell = ctx->cell_start; cell < ctx->cell_end; cell++) {

        manhattan_neighbors_7(cl, cell, neigh);

        const int base_cell = cell * CELL_CAP;

        for (int a = 0; a < cl->counts[cell]; a++) {
            const int i = cl->cells[base_cell + a];
            if (i < 0) continue;

            for (int n = 0; n < 7; n++) {
                const int nb = neigh[n];

                // Cell-level half-shell: only compute each cell-pair once
                if (nb < cell) continue;

                const int base_nb = nb * CELL_CAP;

                for (int b = 0; b < cl->counts[nb]; b++) {
                    const int j = cl->cells[base_nb + b];
                    if (j < 0) continue;

                    // Same-cell half-shell: avoid double counting within cell
                    if (nb == cell && j <= i) continue;

                    double dx = p[i].x - p[j].x;
                    double dy = p[i].y - p[j].y;
                    double dz = p[i].z - p[j].z;

                    md_minimage(&dx, L);
                    md_minimage(&dy, L);
                    md_minimage(&dz, L);

                    const double r2 = dx*dx + dy*dy + dz*dz;
                    if (r2 < tiny2 || r2 > rc2) continue;

                    // Lennard-Jones
                    const double inv_r2  = 1.0 / r2;
                    const double sig2_r2 = (sig * sig) * inv_r2;
                    const double sr6     = sig2_r2 * sig2_r2 * sig2_r2;
                    const double sr12    = sr6 * sr6;

                    U += 4.0 * eps * (sr12 - sr6);

                    const double F_over_r = 24.0 * eps * (2.0 * sr12 - sr6) * inv_r2;

                    const double fx = F_over_r * dx;
                    const double fy = F_over_r * dy;
                    const double fz = F_over_r * dz;

                    // Thread-local Newton updates (no atomics)
                    ctx->local_fx[i] += fx;  ctx->local_fy[i] += fy;  ctx->local_fz[i] += fz;
                    ctx->local_fx[j] -= fx;  ctx->local_fy[j] -= fy;  ctx->local_fz[j] -= fz;
                }
            }
        }
    }

    ctx->local_U = U;
    return NULL;
}

double md_compute_forces_cell_manhattan_pthread(
    Particle *p,
    const SimParams *sp,
    const CellList *cl,
    int nthreads
) {
    if (nthreads < 1) nthreads = 1;
    if (nthreads > cl->ncell) nthreads = cl->ncell;

    const size_t N = sp->N;

    // Zero global forces
    for (size_t i = 0; i < N; i++) {
        p[i].fx = 0.0;
        p[i].fy = 0.0;
        p[i].fz = 0.0;
    }

    pthread_t *threads = (pthread_t*)malloc((size_t)nthreads * sizeof(pthread_t));
    manh_thread_data_t *td = (manh_thread_data_t*)malloc((size_t)nthreads * sizeof(manh_thread_data_t));

    const int cells_per_thread = (cl->ncell + nthreads - 1) / nthreads;

    for (int t = 0; t < nthreads; t++) {
        td[t].cell_start = t * cells_per_thread;
        td[t].cell_end   = td[t].cell_start + cells_per_thread;
        if (td[t].cell_end > cl->ncell) td[t].cell_end = cl->ncell;

        td[t].p  = p;
        td[t].sp = sp;
        td[t].cl = cl;

        td[t].local_fx = (double*)calloc(N, sizeof(double));
        td[t].local_fy = (double*)calloc(N, sizeof(double));
        td[t].local_fz = (double*)calloc(N, sizeof(double));
        td[t].local_U  = 0.0;

        pthread_create(&threads[t], NULL, thread_compute_force_manhattan, &td[t]);
    }

    double U = 0.0;
    for (int t = 0; t < nthreads; t++) {
        pthread_join(threads[t], NULL);
        U += td[t].local_U;
    }

    // Reduce thread-local forces
    for (int t = 0; t < nthreads; t++) {
        for (size_t i = 0; i < N; i++) {
            p[i].fx += td[t].local_fx[i];
            p[i].fy += td[t].local_fy[i];
            p[i].fz += td[t].local_fz[i];
        }
        free(td[t].local_fx);
        free(td[t].local_fy);
        free(td[t].local_fz);
    }

    free(td);
    free(threads);
    return U;
}