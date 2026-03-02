#ifndef CELL_MT_H
#define CELL_MT_H

#include "md.h"
#include "cell_list.h"

// thread local context
struct thread_data {
    int cell_start;
    int cell_end;
    int t_id;
    const CellList *cl;
    Particle *p;
    const SimParams *sp;
    
    // thread-local force buffers
    double *local_fx;
    double *local_fy;
    double *local_fz;
    
    double local_U;
};

double md_compute_forces_cell_mt(Particle *p,
                                const SimParams *sp,
                                const CellList *cl);

#endif