#pragma once

#include <stddef.h>
#include "cell_list.h"
#include "md.h"

// Manhattan (radius-1) cell-list force compute using pthreads.
// - p: particle array (positions + forces)
// - sp: sim params (needs L, rc2, sigma, epsilon, N)
// - cl: built cell list for the CURRENT positions
// - nthreads: number of threads to use (manual)
// Returns potential energy U.
double md_compute_forces_cell_manhattan_pthread(
    Particle *p,
    const SimParams *sp,
    const CellList *cl,
    int nthreads
);