#ifndef PDB_IMPORTER_H
#define PDB_IMPORTER_H

#include "md.h"

//
// Read ATOM lines from a PDB file.
// Returns number of atoms read.
// Also computes bounding box.
//
int pdb_importer(const char *path,
                 Particle *p, int maxN,
                 double *xmin, double *xmax,
                 double *ymin, double *ymax,
                 double *zmin, double *zmax);


int binary_importer(const char *path,
                 Particle *p, int maxN,
                 double *xmin, double *xmax,
                 double *ymin, double *ymax,
                 double *zmin, double *zmax);
#endif
