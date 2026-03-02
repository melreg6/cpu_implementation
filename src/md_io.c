#include "md_io.h"
#include <stdio.h>





void io_write_pdb(const char *path, Particle *p, size_t N)
{
    FILE *f = fopen(path, "w");
    if (!f) {
        perror("fopen");
        return;
    }

    for (size_t i = 0; i < N; i++) {
        fprintf(f,
            "ATOM  %5zu  X   XXX A   1    %8.3f %8.3f %8.3f  1.00  0.00           X\n",
            i + 1,
            p[i].x, p[i].y, p[i].z);
    }

    fprintf(f, "END\n");
    fclose(f);
}
