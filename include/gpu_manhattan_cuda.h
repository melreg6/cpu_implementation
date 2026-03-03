#pragma once
#include "md.h"
#include "cell_list.h"

#ifdef __cplusplus
extern "C" {
#endif

int gpu_manhattan_cuda_init(void);
void gpu_manhattan_cuda_cleanup(void);

/**
 * Compute Lennard-Jones Manhattan cell-list forces on GPU.
 * CPU must have built cl->counts/cl->cells for current positions.
 * Returns potential energy U (computed on GPU + reduced on CPU).
 */
double gpu_manhattan_cuda_forces_and_U(Particle *p, const SimParams *sp, const CellList *cl);

#ifdef __cplusplus
}
#endif