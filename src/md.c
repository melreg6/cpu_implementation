#include "md.h"
#include "cell_list.h"
#include "neighbor_list.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cell_mt.h"

//
// ======================================================
//  FULL O(N^2) FORCE KERNEL
// ======================================================
//
double md_compute_forces_full(Particle *p, const SimParams *sp)
{
    size_t N = sp->N;
    double eps = sp->epsilon;
    double sig = sp->sigma;
    double L   = sp->L;
    double rc2 = sp->rc2;

    const double tiny2 = 1e-30;

    // Reset forces
    for (size_t i = 0; i < N; i++)
        p[i].fx = p[i].fy = p[i].fz = 0.0;

    double U = 0.0;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i + 1; j < N; j++)
        {
            double dx = p[i].x - p[j].x;
            double dy = p[i].y - p[j].y;
            double dz = p[i].z - p[j].z;

            md_minimage(&dx, L);
            md_minimage(&dy, L);
            md_minimage(&dz, L);

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < tiny2 || r2 > rc2)
                continue;

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
// ======================================================
//  FULL O(N^2) VELOCITY VERLET INTEGRATOR
// ======================================================
//
double md_integrate_full(Particle *p, const SimParams *sp, double *Kout)
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
        if (p[i].x < 0) p[i].x += L; else if (p[i].x >= L) p[i].x -= L;
        if (p[i].y < 0) p[i].y += L; else if (p[i].y >= L) p[i].y -= L;
        if (p[i].z < 0) p[i].z += L; else if (p[i].z >= L) p[i].z -= L;
    }

    // New forces
    double U = md_compute_forces_full(p, sp);

    // Final half-step + compute kinetic energy
    double K = 0.0;
    for (size_t i = 0; i < N; i++)
    {
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


//
// ======================================================
//  CELL-LIST MD VELOCITY VERLET (wrapper)
//  Uses forces from md_compute_forces_cell()
// ======================================================
//
double md_integrate_cell(Particle *p, const SimParams *sp,
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
        if (p[i].x < 0) p[i].x += L; else if (p[i].x >= L) p[i].x -= L;
        if (p[i].y < 0) p[i].y += L; else if (p[i].y >= L) p[i].y -= L;
        if (p[i].z < 0) p[i].z += L; else if (p[i].z >= L) p[i].z -= L;
    }

    // Must rebuild cell list every step (cheap)
    cell_list_build(cl, p, L);

    // New forces from cell-list LJ kernel
    double U = md_compute_forces_cell(p, sp, cl);

    // Final half-step + kinetic energy
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

// ======================================================
//  CELL-LIST MD VELOCITY VERLET (HALF-SHELL PARALLEL)
//  Uses forces from md_compute_forces_cell_mt
// ======================================================
double md_integrate_cell_mt(Particle *p, const SimParams *sp,
                                       CellList *cl,
                                       double *Kout)
{
    size_t N = sp->N;
    double dt = sp->dt;
    double half = 0.5 * dt;
    double L = sp->L;

    for (size_t i = 0; i < N; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;
    }

    for (size_t i = 0; i < N; i++) {
        p[i].x += dt * p[i].vx;
        p[i].y += dt * p[i].vy;
        p[i].z += dt * p[i].vz;

        if (p[i].x < 0) p[i].x += L; 
        else if (p[i].x >= L) p[i].x -= L;
        
        if (p[i].y < 0) p[i].y += L; 
        else if (p[i].y >= L) p[i].y -= L;
        
        if (p[i].z < 0) p[i].z += L; 
        else if (p[i].z >= L) p[i].z -= L;
    }

    cell_list_build(cl, p, L);

    // this function handles its own thread creation/joining internally
    double U = md_compute_forces_cell_mt(p, sp, cl);

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

