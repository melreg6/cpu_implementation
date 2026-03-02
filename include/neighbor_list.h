#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H

#include "md.h"
#include "cell_list.h"

//
// Proper named struct
//
typedef struct NeighborList {
    int *nb;
    int *nb_index;
    int total;

    double skin;

    Particle *prev;
} NeighborList;

void nbl_init(NeighborList *nl, size_t N, double rc, double skin);
void nbl_copy_positions_N(NeighborList *nl, Particle *p, size_t N);
int  nbl_needs_rebuild(NeighborList *nl, Particle *p, size_t N, double L);

void nbl_build(NeighborList *nl,
               CellList *cl,
               Particle *p,
               double L,
               double rc,
               size_t N);

#endif
