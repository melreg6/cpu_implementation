#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "md.h"

// GPU Manhattan (radius-1) cell-list force computation.
// - Reads:  p[i].x, p[i].y, p[i].z
// - Writes: p[i].fx, p[i].fy, p[i].fz
// - Returns: total potential energy U (LJ), with double-counting corrected.
double md_compute_forces_cell_manhattan_gpu(Particle *p, const SimParams *sp);

#ifdef __cplusplus
}
#endif