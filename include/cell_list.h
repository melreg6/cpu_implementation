#ifndef CELL_LIST_H
#define CELL_LIST_H

#include "md.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct CellList {
    size_t N;
    int nc;
    int ncell;
    double cell_size;

    int *counts;
    int *cells;
} CellList;

void cell_list_init(CellList *cl, size_t N, double L, double rc);
void cell_list_build(CellList *cl, Particle *p, double L);

int cell_list_collect_neighbors(const CellList *cl,
                                size_t i,
                                const Particle *p,
                                double L,
                                int *buf);

int cl_get_neighbor_cells(const CellList *cl,
                          int cell_idx,
                          int *neighbor_cells);


#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // CELL_LIST_H
