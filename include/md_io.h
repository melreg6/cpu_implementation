#ifndef MD_IO_H
#define MD_IO_H

#include "md.h"

//
// Write particles to PDB file
//
void io_write_pdb(const char *path, Particle *p, size_t N);

#endif
