#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "config.h"
#include "md.h"
#include "md_io.h"
#include "pdb_importer.h"
#include "cell_list.h"
#include "gpu_manhattan_cuda.h"

static double interval(struct timespec a, struct timespec b) {
    struct timespec t;
    t.tv_sec  = b.tv_sec  - a.tv_sec;
    t.tv_nsec = b.tv_nsec - a.tv_nsec;
    if (t.tv_nsec < 0) { t.tv_sec -= 1; t.tv_nsec += 1000000000; }
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

static void copy_particles(Particle *dst, const Particle *src, size_t N) {
    for (size_t i=0;i<N;i++) dst[i] = src[i];
}

static FILE *open_energy_csv(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { perror("fopen"); return NULL; }
    fprintf(f, "step,K,U,E\n");
    return f;
}

static inline double wrap_pos(double x, double L) {
    x = fmod(x, L);
    if (x < 0) x += L;
    return x;
}

static double kinetic_energy(const Particle *p, size_t N) {
    double K = 0.0;
    for (size_t i=0;i<N;i++)
        K += 0.5*(p[i].vx*p[i].vx + p[i].vy*p[i].vy + p[i].vz*p[i].vz);
    return K;
}

static double vv_step_gpu_manhattan(Particle *p, const SimParams *sp, CellList *cl, double *Kout) {
    const size_t N = sp->N;
    const double dt = sp->dt;
    const double L  = sp->L;

    // half step v + position update
    for (size_t i=0;i<N;i++) {
        p[i].vx += 0.5*dt*p[i].fx;
        p[i].vy += 0.5*dt*p[i].fy;
        p[i].vz += 0.5*dt*p[i].fz;

        p[i].x += dt*p[i].vx;
        p[i].y += dt*p[i].vy;
        p[i].z += dt*p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    // rebuild cell list on CPU
    cell_list_build(cl, p, L);

    // GPU forces + potential
    double U = gpu_manhattan_cuda_forces_and_U(p, sp, cl);

    // second half step v
    for (size_t i=0;i<N;i++) {
        p[i].vx += 0.5*dt*p[i].fx;
        p[i].vy += 0.5*dt*p[i].fy;
        p[i].vz += 0.5*dt*p[i].fz;
    }

    if (Kout) *Kout = kinetic_energy(p, N);
    return U;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s input_path\n", argv[0]);
        return 1;
    }
    const char *input_path = argv[1];

    printf("\n================ INITIALIZING GPU MANHATTAN ================\n");

    // Converter (Colab-friendly): python3 (requires pip install MDAnalysis)
    char command[512];
    snprintf(command, sizeof(command),
             "python3 src/py/converter.py %s temp_coords.bin temp_metadata.txt",
             input_path);
    int status = system(command);
    if (status != 0) {
        fprintf(stderr, "Error: Python conversion failed.\n");
        return 1;
    }

    FILE *f_meta = fopen("temp_metadata.txt", "r");
    int N = 0;
    if (!f_meta) { fprintf(stderr, "Metadata not found!\n"); return 1; }
    fscanf(f_meta, "%d", &N);
    fclose(f_meta);

    double xmin, xmax, ymin, ymax, zmin, zmax;
    Particle *p0 = (Particle*)malloc((size_t)N * sizeof(Particle));
    if (!p0) { fprintf(stderr, "malloc p0 failed\n"); return 1; }

    binary_importer("temp_coords.bin", p0, N, &xmin,&xmax,&ymin,&ymax,&zmin,&zmax);

    system("rm -f temp_coords.bin temp_metadata.txt");

    printf("Loaded %d atoms from %s\n", N, input_path);

    // SimParams from config.h
    SimParams sp;
    sp.N = (size_t)N;
    sp.sigma = CONF_SIGMA;
    sp.epsilon = CONF_EPSILON;
    sp.rc = CONF_CUTOFF;
    sp.rc2 = sp.rc * sp.rc;
    sp.dt = CONF_DELTA_T;

    sp.L = fmax(fmax(xmax-xmin, ymax-ymin), zmax-zmin);
    if (sp.L <= 0.0) sp.L = CONF_BOX_MAX;

    printf("Simulation box L = %.3f\n", sp.L);
    printf("rc = %.3f, dt = %.5f\n", sp.rc, sp.dt);

    if (!gpu_manhattan_cuda_init()) {
        fprintf(stderr, "CUDA init failed (no GPU?)\n");
        free(p0);
        return 1;
    }

    // Working copy
    Particle *p = (Particle*)malloc((size_t)N * sizeof(Particle));
    copy_particles(p, p0, (size_t)N);

    // Cell list
    CellList cl;
    cell_list_init(&cl, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl, p, sp.L);

    // initial forces
    double U0 = gpu_manhattan_cuda_forces_and_U(p, &sp, &cl);
    (void)U0;

    // output
    system("mkdir -p output");
    FILE *fE = open_energy_csv("output/gpu_manhattan_energies.csv");
    if (!fE) return 1;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int step=0; step<USER_N_TIMESTEPS; step++) {
        double K, U;
        U = vv_step_gpu_manhattan(p, &sp, &cl, &K);
        fprintf(fE, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double total = interval(t0, t1);

    fclose(fE);
    io_write_pdb("output/gpu_manhattan_positions.pdb", p, (size_t)N);

    printf("\nGPU Manhattan complete.\n");
    printf("Total time: %.6f seconds\n", total);
    printf("Wrote: output/gpu_manhattan_energies.csv\n");
    printf("Wrote: output/gpu_manhattan_positions.pdb\n");

    gpu_manhattan_cuda_cleanup();
    free(cl.counts); free(cl.cells);
    free(p0); free(p);

    return 0;
}