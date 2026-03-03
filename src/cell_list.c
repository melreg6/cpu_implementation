#include "cell_list.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


//
// Helper: convert (ix,iy,iz) to flattened index
//
static inline int cell_index(int ix, int iy, int iz, int nc)
{
    return ix + nc * (iy + nc * iz);
}

//
// Helper: periodic wrap for integer cell coordinates
//
static inline int wrap(int a, int nc)
{
    if (a < 0)     return a + nc;
    if (a >= nc)  return a - nc;
    return a;
}

//
// Initialize cell list structure
// Store N inside the struct (per your choice 1B)
//
void cell_list_init(CellList *cl, size_t N, double L, double rc)
{
    cl->N = N;

    cl->nc = (int)(L / rc);
    if (cl->nc < 1) cl->nc = 1;

    cl->cell_size = L / cl->nc;
    cl->ncell     = cl->nc * cl->nc * cl->nc;

    // counts for each cell
    cl->counts = calloc(cl->ncell, sizeof(int));

    // generous capacity: up to 16 particles per cell (adjustable)
    cl->cells  = calloc(cl->ncell * 16, sizeof(int));
}

//
// Build the cell list for particle array p
//
void cell_list_build(CellList *cl, Particle *p, double L)
{
    int nc = cl->nc;
    size_t N = cl->N;

    // reset per-cell counts
    for (int c = 0; c < cl->ncell; c++)
        cl->counts[c] = 0;

    // clear cell storage
    size_t total_capacity = (size_t)cl->ncell * 16;
    for (size_t idx = 0; idx < total_capacity; idx++)
        cl->cells[idx] = -1;

    double cs = cl->cell_size;

    // assign particles to cells
    for (size_t i = 0; i < N; i++)
    {
        int ix = (int)floor(p[i].x / cs);
        int iy = (int)floor(p[i].y / cs);
        int iz = (int)floor(p[i].z / cs);

        ix = wrap(ix, nc);
        iy = wrap(iy, nc);
        iz = wrap(iz, nc);

        int c = cell_index(ix, iy, iz, nc);

        int k = cl->counts[c]++;
        cl->cells[c * 16 + k] = (int)i;
    }
}

//
// Collect neighbors of particle i by scanning its 27 nearby cells.
// Used by BOTH: cell-list MD + neighbor list builder.
// Caller allocates buf[] of size at least N.
//
int cell_list_collect_neighbors(const CellList *cl,
                                size_t i,
                                const Particle *p,
                                double L,
                                int *buf)
{
    int nc = cl->nc;
    double cs = cl->cell_size;

    // Determine the home cell of particle i
    int ix = (int)floor(p[i].x / cs);
    int iy = (int)floor(p[i].y / cs);
    int iz = (int)floor(p[i].z / cs);

    ix = wrap(ix, nc);
    iy = wrap(iy, nc);
    iz = wrap(iz, nc);

    int count = 0;

    // Loop over the 27 surrounding cells (including home cell)
    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int nx = wrap(ix + dx, nc);
        int ny = wrap(iy + dy, nc);
        int nz = wrap(iz + dz, nc);

        int c = cell_index(nx, ny, nz, nc);
        int base = c * 16;

        for (int k = 0; k < cl->counts[c]; k++) {
            int j = cl->cells[base + k];

            if (j < 0) continue;
            if ((size_t)j == i) continue;

            buf[count++] = j;
        }
    }

    return count;
}

// Collect indicies of neighbor cells for a given cell (including itself)
int cl_get_neighbor_cells(const CellList *cl, int cell_idx, int *neighbor_cells)
{
    int nc = cl->nc;

    // convert 1D index to 3D coordinates
    int z = cell_idx / (nc * nc);
    int rem = cell_idx % (nc * nc);
    int y = rem / nc;
    int x = rem % nc;

    int count = 0;

    // loop over the 3x3x3 neighborhood
    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int nx = wrap(x + dx, nc);
        int ny = wrap(y + dy, nc);
        int nz = wrap(z + dz, nc);

        int neighbor_idx = cell_index(nx, ny, nz, nc);
        neighbor_cells[count++] = neighbor_idx;
    }

    return count;
}

//
// Full MD force computation using cell lists
// (per your choice 2C)
//
double md_compute_forces_cell(Particle *p, const SimParams *sp,
                              const CellList *cl)
{
    size_t N = sp->N;
    double eps = sp->epsilon;
    double sig = sp->sigma;
    double L   = sp->L;
    double rc2 = sp->rc2;

    // zero forces
    for (size_t i = 0; i < N; i++)
        p[i].fx = p[i].fy = p[i].fz = 0.0;

    double U = 0.0;
    const double tiny2 = 1e-30;

    // temp buffer to store neighbors
    int *buf = malloc(N * sizeof(int));

    for (size_t i = 0; i < N; i++)
    {
        int nnb = cell_list_collect_neighbors(cl, i, p, L, buf);

        for (int t = 0; t < nnb; t++) {
            int j = buf[t];
            if (j < 0) continue;

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

            // Newton's third law
            p[j].fx -= fx;
            p[j].fy -= fy;
            p[j].fz -= fz;

            U += 4.0 * eps * (sr12 - sr6);
        }
    }

    free(buf);
    return U;
}
