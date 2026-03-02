#include "neighbor_list.h"
#include "cell_list.h"
#include "md.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//
// Allocate and initialize neighbor list
//
void nbl_init(NeighborList *nl, size_t N, double rc, double skin)
{
    nl->nb       = calloc(N * 64, sizeof(int));   // generous capacity
    nl->nb_index = calloc(N + 1, sizeof(int));
    nl->prev     = calloc(N, sizeof(Particle));

    nl->total    = 0;
    nl->skin     = skin;
}

//
// Copy current particle positions so we know when to rebuild
//
void nbl_copy_positions(NeighborList *nl, Particle *p)
{
    // Only copy x,y,z — no need for velocities
    for (size_t i = 0; i < nl->nb_index[0]; i++) {
        nl->prev[i].x = p[i].x;
        nl->prev[i].y = p[i].y;
        nl->prev[i].z = p[i].z;
    }
}

// Actually: nb_index[0] doesn't store N.
// So we fix that by passing N in from main.

void nbl_copy_positions_N(NeighborList *nl, Particle *p, size_t N)
{
    for (size_t i = 0; i < N; i++) {
        nl->prev[i].x = p[i].x;
        nl->prev[i].y = p[i].y;
        nl->prev[i].z = p[i].z;
    }
}

//
// Check if we need to rebuild neighbor list (has any particle moved > skin/2?)
//
int nbl_needs_rebuild(NeighborList *nl, Particle *p, size_t N, double L)
{
    double thresh = 0.5 * nl->skin;
    double thresh2 = thresh * thresh;

    for (size_t i = 0; i < N; i++) {
        double dx = p[i].x - nl->prev[i].x;
        double dy = p[i].y - nl->prev[i].y;
        double dz = p[i].z - nl->prev[i].z;

        md_minimage(&dx, L);
        md_minimage(&dy, L);
        md_minimage(&dz, L);

        double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 > thresh2)
            return 1;   // REBUILD
    }
    return 0; // no rebuild needed
}

//
// Build neighbor list using the cell list for candidate pairs
//
void nbl_build(NeighborList *nl, CellList *cl, Particle *p, double L, double rc, size_t N)
{
    double rc_skin = rc + nl->skin;
    double rc_skin2 = rc_skin * rc_skin;

    // Reset neighbor list
    nl->total = 0;

    // Reset nb_index
    for (size_t i = 0; i <= N; i++)
        nl->nb_index[i] = 0;

    int *buf = malloc(N * sizeof(int));

    for (size_t i = 0; i < N; i++)
    {
        // Where this particle's neighbor list begins
        nl->nb_index[i] = nl->total;

        // Collect neighbor candidates from 27 cells
        int nnb = cell_list_collect_neighbors(cl, i, p, L, buf);

        for (int t = 0; t < nnb; t++) {
            int j = buf[t];

            double dx = p[i].x - p[j].x;
            double dy = p[i].y - p[j].y;
            double dz = p[i].z - p[j].z;

            md_minimage(&dx, L);
            md_minimage(&dy, L);
            md_minimage(&dz, L);

            double r2 = dx*dx + dy*dy + dz*dz;

            // Include if within (rc + skin)^2
            if (r2 <= rc_skin2) {
                nl->nb[nl->total++] = j;
            }
        }
    }

    nl->nb_index[N] = nl->total;

    free(buf);

    // Finally store current positions
    nbl_copy_positions_N(nl, p, N);

//    printf("[NBL] Built neighbor list: total entries %d (avg %.2f per atom)\n",
 //          nl->total, (double)nl->total / (double)N);
}

//
// Force calculation using neighbor list
//
double md_compute_forces_nbl(Particle *p,
                             const SimParams *sp,
                             const NeighborList *nl)
{
    size_t N = sp->N;
    double eps = sp->epsilon;
    double sig = sp->sigma;
    double L   = sp->L;
    double rc2 = sp->rc2;

    // Zero forces
    for (size_t i = 0; i < N; i++)
        p[i].fx = p[i].fy = p[i].fz = 0.0;

    double U = 0.0;
    const double tiny2 = 1e-30;

    for (size_t i = 0; i < N; i++)
    {
        int start = nl->nb_index[i];
        int end   = nl->nb_index[i+1];

        for (int k = start; k < end; k++) {
            int j = nl->nb[k];

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
            double sr6  = sig2_r2 * sig2_r2 * sig2_r2;
            double sr12 = sr6 * sr6;

            double F_over_r = 24.0 * eps * (2.0*sr12 - sr6) * inv_r2;

            double fx = F_over_r * dx;
            double fy = F_over_r * dy;
            double fz = F_over_r * dz;

            p[i].fx += fx;
            p[i].fy += fy;
            p[i].fz += fz;

            p[j].fx -= fx;
            p[j].fy -= fy;
            p[j].fz -= fz;

            U += 4.0 * eps * (sr12 - sr6);
        }
    }

    return U;
}

//
// Velocity Verlet integrator for neighbor-list MD
//
double md_integrate_nbl(Particle *p,
                        const SimParams *sp,
                        NeighborList *nl,
                        CellList *cl,
                        double *Kout)
{
    size_t N = sp->N;
    double dt = sp->dt;
    double half = 0.5 * dt;
    double L = sp->L;

    // Half-step velocities
    for (size_t i = 0; i < N; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;
    }

    // Update positions
    for (size_t i = 0; i < N; i++) {
        p[i].x += dt * p[i].vx;
        p[i].y += dt * p[i].vy;
        p[i].z += dt * p[i].vz;

        // periodic wrap
        if (p[i].x < 0) p[i].x += L;
        if (p[i].x >= L) p[i].x -= L;
        if (p[i].y < 0) p[i].y += L;
        if (p[i].y >= L) p[i].y -= L;
        if (p[i].z < 0) p[i].z += L;
        if (p[i].z >= L) p[i].z -= L;
    }

    // Need rebuild?
    if (nbl_needs_rebuild(nl, p, N, L)) {
//        printf("[NBL] Rebuilding...\n");
        cell_list_build(cl, p, L);
        nbl_build(nl, cl, p, L, sp->rc, N);
    }

    // Compute new forces
    double U = md_compute_forces_nbl(p, sp, nl);

    // Final half-step velocities + kinetic energy
    double K = 0.0;
    for (size_t i = 0; i < N; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;

        K += 0.5 * (p[i].vx*p[i].vx +
                    p[i].vy*p[i].vy +
                    p[i].vz*p[i].vz);
    }

    *Kout = K;
    return U;
}
