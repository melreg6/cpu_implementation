#ifndef NEIGHBOR_LIST_PTHREAD_H
#define NEIGHBOR_LIST_PTHREAD_H

#include "md.h"
#include "cell_list.h"
#include "neighbor_list.h"
#include <stddef.h>

/**
 * Half-Shell Neighbor List with Pthreads - LARGE SYSTEM OPTIMIZED
 * 
 * Compile with: -pthread -DNBL_NUM_THREADS=8
 */

void nbl_pthread_init(size_t N);
void nbl_pthread_cleanup(void);

void nbl_build_pthread(NeighborList *nl, CellList *cl, Particle *p,
                       double L, double rc, size_t N);

double md_compute_forces_nbl_pthread(Particle *p,
                                     const SimParams *sp,
                                     const NeighborList *nl);

double md_integrate_nbl_pthread(Particle *p,
                                const SimParams *sp,
                                NeighborList *nl,
                                CellList *cl,
                                double *Kout);

int nbl_pthread_get_num_threads(void);
void nbl_pthread_print_info(void);
size_t nbl_pthread_memory_usage(size_t N);

#endif /* NEIGHBOR_LIST_PTHREAD_H */