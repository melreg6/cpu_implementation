#include <stdlib.h>
#include <pthread.h>
#include <string.h> 
#include "cell_mt.h"
#include "cell_list.h"
#include "md.h"

/*
 * thread routine using half shell logic.
 * each thread writes to its own private local_fx/y/z arrays.
 */
static void* thread_compute_force_cell(void *arg) {
    struct thread_data *ctx = (struct thread_data*)arg;
    const CellList *cl = ctx->cl;
    Particle *p = ctx->p;
    const SimParams *sp = ctx->sp;

    double rc2 = sp->rc2;
    double eps = sp->epsilon;
    double sig = sp->sigma;
    double L   = sp->L;
    const double tiny2 = 1e-30;

    double local_U = 0.0;
    int max_cell_particles = 16;
    int neighbor_cells[27];

    // iterate over cells assigned to this thread
    for (int cell = ctx->cell_start; cell < ctx->cell_end; cell++) {
        int base_cell = cell * max_cell_particles;

        // get all 27 neighbor cell indicies (including itself)
        int n_neighbors = cl_get_neighbor_cells(cl, cell, neighbor_cells);

        // iterate over particles in the current cell
        for (int a = 0; a < cl->counts[cell]; a++) {
            int i = cl->cells[base_cell + a];
            if (i < 0) continue;

            // loop over neighbor cells
            for (int n = 0; n < n_neighbors; n++) {
                int neighbor_cell = neighbor_cells[n];

                // cell-level half-shell condition
                if (neighbor_cell < cell) continue;  // skip lower-index cells

                int base_neighbor = neighbor_cell * max_cell_particles;

                // loop over particles in neighbor cell
                for (int b = 0; b < cl->counts[neighbor_cell]; b++) {
                    int j = cl->cells[base_neighbor + b];
                    if (j < 0) continue;

                    // particle-level half-shell for self-cell
                    if (neighbor_cell == cell && j <= i) continue;

                    // compute displacement
                    double dx = p[i].x - p[j].x;
                    double dy = p[i].y - p[j].y;
                    double dz = p[i].z - p[j].z;

                    md_minimage(&dx, L);
                    md_minimage(&dy, L);
                    md_minimage(&dz, L);

                    double r2 = dx*dx + dy*dy + dz*dz;
                    if (r2 < tiny2 || r2 > rc2) continue;

                    double inv_r2 = 1.0 / r2;
                    double sig2_r2 = (sig*sig) * inv_r2;
                    double sr6 = sig2_r2 * sig2_r2 * sig2_r2;
                    double sr12 = sr6 * sr6;
                    double F_over_r = 24.0 * eps * (2.0*sr12 - sr6) * inv_r2;

                    double fx = F_over_r * dx;
                    double fy = F_over_r * dy;
                    double fz = F_over_r * dz;

                    // update thread-local forces
                    ctx->local_fx[i] += fx;
                    ctx->local_fy[i] += fy;
                    ctx->local_fz[i] += fz;

                    ctx->local_fx[j] -= fx;
                    ctx->local_fy[j] -= fy;
                    ctx->local_fz[j] -= fz;

                    local_U += 4.0 * eps * (sr12 - sr6);
                }
            }
        }
    }
    ctx->local_U = local_U;
    return NULL;
}


double md_compute_forces_cell_mt(Particle *p, const SimParams *sp, const CellList *cl) {
    int nthreads = sp->nthreads;
    pthread_t *threads = malloc(nthreads * sizeof(pthread_t));
    struct thread_data *td = malloc(nthreads * sizeof(struct thread_data));

    // init thread data and private buffers
    int cells_per_thread = (cl->ncell + nthreads - 1) / nthreads;
    
    for (int t = 0; t < nthreads; t++) {
        td[t].t_id = t;
        td[t].cell_start = t * cells_per_thread;
        td[t].cell_end = (td[t].cell_start + cells_per_thread > cl->ncell) ? cl->ncell : td[t].cell_start + cells_per_thread;
        
        td[t].cl = cl;
        td[t].p = p;
        td[t].sp = sp;
        td[t].local_U = 0.0;

        td[t].local_fx = calloc(sp->N, sizeof(double));
        td[t].local_fy = calloc(sp->N, sizeof(double));
        td[t].local_fz = calloc(sp->N, sizeof(double));

        pthread_create(&threads[t], NULL, thread_compute_force_cell, &td[t]);
    }

    // wait for computation to finish
    double total_U = 0.0;
    for (int t = 0; t < nthreads; t++) {
        pthread_join(threads[t], NULL);
        total_U += td[t].local_U;
    }

    for (size_t i = 0; i < sp->N; i++) {
        p[i].fx = p[i].fy = p[i].fz = 0.0;
    }

    // sum thread-local forces into global particle forces
    for (int t = 0; t < nthreads; t++) {
        for (size_t i = 0; i < sp->N; i++) {
            p[i].fx += td[t].local_fx[i];
            p[i].fy += td[t].local_fy[i];
            p[i].fz += td[t].local_fz[i];
        }
    }

    // clean
    for (int t = 0; t < nthreads; t++) {
        free(td[t].local_fx);
        free(td[t].local_fy);
        free(td[t].local_fz);
    }
    free(td);
    free(threads);

    return total_U;
}
