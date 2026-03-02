#ifndef MD_H
#define MD_H

#include <stddef.h>

//
// Particle structure — OLD FORMAT
//
typedef struct {
    double x, y, z;        
    double vx, vy, vz;     
    double fx, fy, fz;     
} Particle;

//
// Simulation parameters
//
typedef struct {
    size_t N;
    double L;

    double sigma;
    double epsilon;

    double rc;
    double rc2;

    double dt;

    int nthreads;
} SimParams;

//
// Minimum-image periodic BC
//
static inline void md_minimage(double *d, double L)
{
    if (*d >  0.5*L) *d -= L;
    if (*d < -0.5*L) *d += L;
}

//
// Forward declarations ONLY (no typedefs here)
//
struct CellList;
struct NeighborList;

//
// FULL MD
//
double md_compute_forces_full(Particle *p, const SimParams *sp);
double md_integrate_full(Particle *p, const SimParams *sp, double *Kout);

//
// CELL-LIST MD
//
double md_compute_forces_cell(Particle *p, const SimParams *sp,
                              const struct CellList *cl);

double md_integrate_cell(Particle *p, const SimParams *sp,
                         struct CellList *cl,
                         double *Kout);

//
// NEIGHBOR-LIST MD
//
double md_compute_forces_nbl(Particle *p, const SimParams *sp,
                             const struct NeighborList *nl);

double md_integrate_nbl(Particle *p, const SimParams *sp,
                        struct NeighborList *nl,
                        struct CellList *cl,
                        double *Kout);

//
// CELL LIST MULTITHREADED MD
//
double md_integrate_cell_mt(Particle *p, const SimParams *sp,
                                       struct CellList *cl,
                                       double *Kout);

#endif
