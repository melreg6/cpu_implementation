#ifndef CELL_LIST_H
#define CELL_LIST_H

#include "md.h"

//
// Proper named struct definition (fixes typedef conflict)
//
typedef struct CellList {
    size_t N;        
    int nc;          
    int ncell;       
    double cell_size;

    int *counts;
    int *cells;
} CellList;

//
// API
//
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

double md_compute_forces_cell(Particle *p, const SimParams *sp,
                              const CellList *cl);

#endif
