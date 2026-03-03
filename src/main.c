#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>  // For clock_gettime

#include "config.h"
#include "md.h"
#include "md_io.h"
#include "pdb_importer.h"
#include "cell_list.h"
#include "neighbor_list.h"
#include "cell_mt.h"
#include "neighbor_list_pthread.h"
#include "cell_manhattan_pthread.h"
#include "gpu_manhattan.h"

#define MAX_ATOMS 50000

// ----------------------------------------------------
// Time measurement by clock_gettime()
// ----------------------------------------------------
static double interval(struct timespec a, struct timespec b)
{
    return (b.tv_sec - a.tv_sec) + 1e-9*(b.tv_nsec - a.tv_nsec);
}

// ----------------------------------------------------
// Helpers
// ----------------------------------------------------
static double wrap_pos(double x, double L)
{
    x = fmod(x, L);
    if (x < 0) x += L;
    return x;
}

// Initialize velocities to random in [-0.5, 0.5]
static void init_velocities(Particle *p, size_t N)
{
    for (size_t i = 0; i < N; i++) {
        p[i].vx = ((double)rand() / RAND_MAX) - 0.5;
        p[i].vy = ((double)rand() / RAND_MAX) - 0.5;
        p[i].vz = ((double)rand() / RAND_MAX) - 0.5;
    }
}

static FILE* open_energy_csv(const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) {
        perror("fopen");
        return NULL;
    }
    fprintf(f, "step,K,U,E\n");
    return f;
}

// ----------------------------------------------------
// Integrators (Velocity Verlet style)
// ----------------------------------------------------
static double md_integrate_full(
    Particle *p,
    const SimParams *sp,
    double *K_out
){
    const size_t N = sp->N;
    const double dt = sp->dt;
    const double L  = sp->L;

    // half-step velocity + position
    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;

        p[i].x  += dt * p[i].vx;
        p[i].y  += dt * p[i].vy;
        p[i].z  += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    // forces
    double U = md_compute_forces_full(p, sp);

    // half-step velocity
    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;
    }

    if (K_out) *K_out = kinetic_energy(p, N);
    return U;
}

static double md_integrate_cell(
    Particle *p,
    const SimParams *sp,
    CellList *cl,
    double *K_out
){
    const size_t N = sp->N;
    const double dt = sp->dt;
    const double L  = sp->L;

    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;

        p[i].x  += dt * p[i].vx;
        p[i].y  += dt * p[i].vy;
        p[i].z  += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    cell_list_build(cl, p, L);

    double U = md_compute_forces_cell(p, sp, cl);

    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;
    }

    if (K_out) *K_out = kinetic_energy(p, N);
    return U;
}

static double md_integrate_nbl(
    Particle *p,
    const SimParams *sp,
    NeighborList *nl,
    CellList *cl,
    double *K_out
){
    const size_t N = sp->N;
    const double dt = sp->dt;
    const double L  = sp->L;

    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;

        p[i].x  += dt * p[i].vx;
        p[i].y  += dt * p[i].vy;
        p[i].z  += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    // rebuild bins + neighbor list
    cell_list_build(cl, p, L);
    nbl_build(nl, cl, p, sp->L, sp->rc, sp->N);

    double U = md_compute_forces_nbl(p, sp, nl);

    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;
    }

    if (K_out) *K_out = kinetic_energy(p, N);
    return U;
}

static double md_integrate_cell_mt(
    Particle *p,
    const SimParams *sp,
    CellList *cl,
    double *K_out
){
    const size_t N = sp->N;
    const double dt = sp->dt;
    const double L  = sp->L;

    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;

        p[i].x  += dt * p[i].vx;
        p[i].y  += dt * p[i].vy;
        p[i].z  += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    cell_list_build(cl, p, L);

    double U = md_compute_forces_cell_mt(p, sp, cl);

    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;
    }

    if (K_out) *K_out = kinetic_energy(p, N);
    return U;
}

static double md_integrate_nbl_pthread(
    Particle *p,
    const SimParams *sp,
    NeighborList *nl,
    CellList *cl,
    double *K_out
){
    const size_t N = sp->N;
    const double dt = sp->dt;
    const double L  = sp->L;

    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;

        p[i].x  += dt * p[i].vx;
        p[i].y  += dt * p[i].vy;
        p[i].z  += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    cell_list_build(cl, p, L);
    nbl_build_pthread(nl, cl, p, sp->L, sp->rc, sp->N);

    double U = md_compute_forces_nbl_pthread(p, sp, nl);

    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;
    }

    if (K_out) *K_out = kinetic_energy(p, N);
    return U;
}

static double md_integrate_manhattan_pthread(
    Particle *p,
    const SimParams *sp,
    CellList *cl,
    double *K_out,
    int nthreads
){
    const size_t N = sp->N;
    const double dt = sp->dt;
    const double L  = sp->L;

    // 1) half-step velocity update + position update
    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;

        p[i].x  += dt * p[i].vx;
        p[i].y  += dt * p[i].vy;
        p[i].z  += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    // 2) rebuild bins for new positions
    cell_list_build(cl, p, L);

    // 3) compute new forces + potential
    double U = md_compute_forces_cell_manhattan_pthread(p, sp, cl, nthreads);

    // 4) half-step velocity update using new forces
    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;
    }

    // 5) kinetic energy
    if (K_out) *K_out = kinetic_energy(p, N);
    return U;
}

static double md_integrate_manhattan_gpu(
    Particle *p,
    const SimParams *sp,
    double *K_out
){
    const size_t N = sp->N;
    const double dt = sp->dt;
    const double L  = sp->L;

    // 1) half-step velocity update + position update
    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;

        p[i].x  += dt * p[i].vx;
        p[i].y  += dt * p[i].vy;
        p[i].z  += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    // 2) compute new forces + potential on GPU (GPU builds Manhattan cell bins internally)
    double U = md_compute_forces_cell_manhattan_gpu(p, sp);

    // 3) half-step velocity update using new forces
    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;
    }

    // 4) kinetic energy
    if (K_out) *K_out = kinetic_energy(p, N);
    return U;
}

// ---------------------------------------
// MAIN
// ---------------------------------------
int main(int argc, char **argv)
{
    srand(0);

    // ------------------------------------------------
    // Parse args
    // ------------------------------------------------
    if (argc < 2) {
        printf("Usage: %s input.pdb [--threads N]\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    int user_threads = 4;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            user_threads = atoi(argv[i + 1]);
            i++;
        }
    }

    // ------------------------------------------------
    // Read PDB
    // ------------------------------------------------
    Particle p0[MAX_ATOMS];
    int N = 0;
    if (pdb_import_particles(input_path, p0, &N) != 0) {
        fprintf(stderr, "Failed to read input PDB: %s\n", input_path);
        return 1;
    }
    printf("Read %d particles from %s\n", N, input_path);

    // ------------------------------------------------
    // Params
    // ------------------------------------------------
    SimParams sp;
    sp.N = (size_t)N;
    sp.L = 50.0;           // box length (tune as needed)
    sp.rc = 2.5;
    sp.rc2 = sp.rc * sp.rc;
    sp.sigma = 1.0;
    sp.epsilon = 1.0;
    sp.dt = 0.005;
    sp.nthreads = user_threads;

    printf("Box L=%.2f, rc=%.2f, dt=%.4f\n", sp.L, sp.rc, sp.dt);

    // ------------------------------------------------
    // Make copies for each algorithm
    // ------------------------------------------------
    Particle *p0_dyn = (Particle*)malloc((size_t)N * sizeof(Particle));
    Particle *p_full = (Particle*)malloc((size_t)N * sizeof(Particle));
    Particle *p_cell = (Particle*)malloc((size_t)N * sizeof(Particle));
    Particle *p_nbl  = (Particle*)malloc((size_t)N * sizeof(Particle));
    Particle *p_cell_mt = (Particle*)malloc((size_t)N * sizeof(Particle));
    Particle *p_nbl_pthread = (Particle*)malloc((size_t)N * sizeof(Particle));
    Particle *p_manhattan = (Particle*)malloc((size_t)N * sizeof(Particle));
    Particle *p_final = (Particle*)malloc((size_t)N * sizeof(Particle));

    if (!p0_dyn || !p_full || !p_cell || !p_nbl || !p_cell_mt || !p_nbl_pthread || !p_manhattan || !p_final) {
        fprintf(stderr, "malloc failed\n");
        free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        p0_dyn[i] = p0[i];
    }

    // initialize velocities + forces
    init_velocities(p0_dyn, (size_t)N);
    zero_forces(p0_dyn, (size_t)N);

    // copies
    memcpy(p_full, p0_dyn, (size_t)N * sizeof(Particle));
    memcpy(p_cell, p0_dyn, (size_t)N * sizeof(Particle));
    memcpy(p_nbl,  p0_dyn, (size_t)N * sizeof(Particle));
    memcpy(p_cell_mt, p0_dyn, (size_t)N * sizeof(Particle));
    memcpy(p_nbl_pthread, p0_dyn, (size_t)N * sizeof(Particle));
    memcpy(p_manhattan, p0_dyn, (size_t)N * sizeof(Particle));
    memcpy(p_final, p0_dyn, (size_t)N * sizeof(Particle));

    // ------------------------------------------------
    // Setup output dir
    // ------------------------------------------------
    system("mkdir -p output");

    // ------------------------------------------------
    // Init common data structures
    // ------------------------------------------------
    CellList cl;
    cell_list_init(&cl, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl, p_cell, sp.L);

    CellList cl2;
    cell_list_init(&cl2, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl2, p_nbl, sp.L);

    NeighborList nl;
    nbl_init(&nl, (size_t)N, sp.rc, 0.3 * sp.rc);
    nbl_build(&nl, &cl2, p_nbl, sp.L, sp.rc, (size_t)N);

    // pthread neighbor list system init
    nbl_pthread_init(sp.nthreads);

    // ------------------------------------------------
    // AUTOTUNE SETTINGS
    // ------------------------------------------------
   const int AUTOTUNE_N_TIMESTEPS = 50;

    struct timespec time_start, time_stop;

    // ------------------------------------------------
    // 1) FULL O(N^2)
    // ------------------------------------------------
    printf("\n================ FULL O(N^2) ================\n");
    FILE *f_full = open_energy_csv("output/full_energies.csv");
    if (!f_full) {
        free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    md_compute_forces_full(p_full, &sp);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_full(p_full, &sp, &K);
        fprintf(f_full, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_full = interval(time_start, time_stop);

    fclose(f_full);
    io_write_pdb("output/full_positions.pdb", p_full, (size_t)N);
    printf("FULL done. Time: %.6f seconds\n", time_full);

    // ------------------------------------------------
    // 2) CELL-LIST
    // ------------------------------------------------
    printf("\n================ CELL-LIST ================\n");
    FILE *f_cell = open_energy_csv("output/cell_energies.csv");
    if (!f_cell) {
        free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    md_compute_forces_cell(p_cell, &sp, &cl);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_cell(p_cell, &sp, &cl, &K);
        fprintf(f_cell, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_cell = interval(time_start, time_stop);

    fclose(f_cell);
    io_write_pdb("output/cell_positions.pdb", p_cell, (size_t)N);
    printf("CELL-LIST done. Time: %.6f seconds\n", time_cell);

    // ------------------------------------------------
    // 3) NEIGHBOR-LIST (Serial)
    // ------------------------------------------------
    printf("\n================ NEIGHBOR-LIST (Serial) ================\n");
    FILE *f_nbl = open_energy_csv("output/nbl_energies.csv");
    if (!f_nbl) {
        free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    md_compute_forces_nbl(p_nbl, &sp, &nl);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_nbl(p_nbl, &sp, &nl, &cl2, &K);
        fprintf(f_nbl, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_nbl = interval(time_start, time_stop);

    fclose(f_nbl);
    io_write_pdb("output/nbl_positions.pdb", p_nbl, (size_t)N);
    printf("NEIGHBOR-LIST (Serial) done. Time: %.6f seconds\n", time_nbl);

    // ------------------------------------------------
    // 4) CELL-LIST HALF-SHELL MT
    // ------------------------------------------------
    printf("\n================ CELL-LIST HALF-SHELL MT ================\n");
    FILE *f_cell_mt = open_energy_csv("output/cell_mt_energies.csv");
    if (!f_cell_mt) {
        free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    CellList cl_mt;
    cell_list_init(&cl_mt, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl_mt, p_cell_mt, sp.L);

    md_compute_forces_cell_mt(p_cell_mt, &sp, &cl_mt);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_cell_mt(p_cell_mt, &sp, &cl_mt, &K);
        fprintf(f_cell_mt, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_cell_mt = interval(time_start, time_stop);

    fclose(f_cell_mt);
    io_write_pdb("output/cell_mt_positions.pdb", p_cell_mt, (size_t)N);
    printf("CELL-LIST HALF-SHELL MT done. Time: %.6f seconds\n", time_cell_mt);

    // ------------------------------------------------
    // 5) NEIGHBOR-LIST (Pthread)
    // ------------------------------------------------
    printf("\n================ NEIGHBOR-LIST (Pthread) ================\n");
    printf("NBL pthread threads: %d\n", nbl_pthread_get_num_threads());

    FILE *f_nbl_pthread = open_energy_csv("output/nbl_pthread_energies.csv");
    if (!f_nbl_pthread) {
        free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    CellList cl_nbl_pt;
    cell_list_init(&cl_nbl_pt, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl_nbl_pt, p_nbl_pthread, sp.L);

    NeighborList nl_pt;
    nbl_init(&nl_pt, (size_t)N, sp.rc, 0.3 * sp.rc);
    nbl_build_pthread(&nl_pt, &cl_nbl_pt, p_nbl_pthread, sp.L, sp.rc, (size_t)N);

    md_compute_forces_nbl_pthread(p_nbl_pthread, &sp, &nl_pt);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_nbl_pthread(p_nbl_pthread, &sp, &nl_pt, &cl_nbl_pt, &K);
        fprintf(f_nbl_pthread, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_nbl_pthread = interval(time_start, time_stop);

    fclose(f_nbl_pthread);
    io_write_pdb("output/nbl_pthread_positions.pdb", p_nbl_pthread, (size_t)N);
    printf("NEIGHBOR-LIST (Pthread) done. Time: %.6f seconds\n", time_nbl_pthread);

    // ------------------------------------------------
    // 7.5 MANHATTAN CELL-LIST MD (Pthread)
    // ------------------------------------------------
    printf("\n================ MANHATTAN CELL-LIST MD (Pthread) ================\n");
    printf("Manhattan threads: %d\n", sp.nthreads);

    FILE *f_manh = open_energy_csv("output/manhattan_pthread_energies.csv");
    if (!f_manh) {
        free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    CellList cl_manh;
    cell_list_init(&cl_manh, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl_manh, p_manhattan, sp.L);

    // initial forces
    md_compute_forces_cell_manhattan_pthread(p_manhattan, &sp, &cl_manh, sp.nthreads);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_manhattan_pthread(p_manhattan, &sp, &cl_manh, &K, sp.nthreads);
        fprintf(f_manh, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_manh = interval(time_start, time_stop);

    fclose(f_manh);
    io_write_pdb("output/manhattan_pthread_positions.pdb", p_manhattan, (size_t)N);
    printf("MANHATTAN CELL-LIST (Pthread) done. Time: %.6f seconds\n", time_manh);

    // ------------------------------------------------
    // 7.6 MANHATTAN CELL-LIST MD (GPU)
    // ------------------------------------------------
    printf("\n================ MANHATTAN CELL-LIST MD (GPU) ================\n");

    FILE *f_manh_gpu = open_energy_csv("output/manhattan_gpu_energies.csv");
    if (!f_manh_gpu) {
        free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    // initial forces on GPU
    md_compute_forces_cell_manhattan_gpu(p_manhattan, &sp);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_manhattan_gpu(p_manhattan, &sp, &K);
        fprintf(f_manh_gpu, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_manh_gpu = interval(time_start, time_stop);

    fclose(f_manh_gpu);
    io_write_pdb("output/manhattan_gpu_positions.pdb", p_manhattan, (size_t)N);
    printf("MANHATTAN CELL-LIST (GPU) done. Time: %.6f seconds\n", time_manh_gpu);

    // ------------------------------------------------
    // 8. COMPARE TIMES AND DETERMINE FASTEST
    // ------------------------------------------------
    printf("\n========================================\n");
    printf("PERFORMANCE COMPARISON\n");
    printf("========================================\n");
    printf("FULL O(N^2):                %.6f seconds (1.00x speedup)\n", time_full);
    printf("CELL-LIST:                  %.6f seconds (%.2fx speedup)\n", time_cell, time_full/time_cell);
    printf("NEIGHBOR-LIST (Serial):     %.6f seconds (%.2fx speedup)\n", time_nbl, time_full/time_nbl);
    printf("CELL-LIST HALF-SHELL MT:    %.6f seconds (%.2fx speedup)\n", time_cell_mt, time_full/time_cell_mt);
    printf("NBL-PTHREAD (%d threads):   %.6f seconds (%.2fx speedup)\n",
           nbl_pthread_get_num_threads(), time_nbl_pthread, time_full/time_nbl_pthread);
    printf("MANHATTAN (%d threads):     %.6f seconds (%.2fx speedup)\n",
           sp.nthreads, time_manh, time_full/time_manh);
    printf("MANHATTAN (GPU):            %.6f seconds (%.2fx speedup)\n",
           time_manh_gpu, time_full/time_manh_gpu);
    printf("========================================\n");

    // Determine fastest
    const char *fastest_name = "FULL O(N^2)";
    int fastest_method = 0;
    double fastest_time = time_full;

    if (time_cell < fastest_time) {
        fastest_time = time_cell;
        fastest_method = 1;
        fastest_name = "CELL-LIST";
    }
    if (time_nbl < fastest_time) {
        fastest_time = time_nbl;
        fastest_method = 2;
        fastest_name = "NEIGHBOR-LIST (Serial)";
    }
    if (time_cell_mt < fastest_time) {
        fastest_time = time_cell_mt;
        fastest_method = 3;
        fastest_name = "CELL-LIST HALF-SHELL MT";
    }
    if (time_nbl_pthread < fastest_time) {
        fastest_time = time_nbl_pthread;
        fastest_method = 4;
        fastest_name = "NBL-PTHREAD";
    }
    if (time_manh < fastest_time) {
        fastest_time = time_manh;
        fastest_method = 5;
        fastest_name = "MANHATTAN CELL-LIST (Pthread)";
    }

    if (time_manh_gpu < fastest_time) {
        fastest_time = time_manh_gpu;
        fastest_method = 6;
        fastest_name = "MANHATTAN CELL-LIST (GPU)";
    }

    printf("\nFASTEST METHOD: %s (%.6f seconds)\n", fastest_name, fastest_time);

    // ------------------------------------------------
    // 9. RUN FULL SIMULATION WITH FASTEST METHOD
    // ------------------------------------------------
    const int USER_N_TIMESTEPS = 200;

    printf("\n========================================\n");
    printf("RUNNING FULL SIMULATION WITH: %s\nfor %d TIMESTEPS\n",
           fastest_name, USER_N_TIMESTEPS);
    printf("========================================\n");

    FILE *f_final = open_energy_csv("output/final_energies.csv");
    if (!f_final) {
        free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    if (fastest_method == 0) {
        md_compute_forces_full(p_final, &sp);
        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_full(p_final, &sp, &K);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }
    }
    else if (fastest_method == 1) {
        CellList cl_final;
        cell_list_init(&cl_final, (size_t)N, sp.L, sp.rc);
        cell_list_build(&cl_final, p_final, sp.L);

        md_compute_forces_cell(p_final, &sp, &cl_final);

        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_cell(p_final, &sp, &cl_final, &K);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }

        free(cl_final.counts);
        free(cl_final.cells);
    }
    else if (fastest_method == 2) {
        CellList cl_final;
        cell_list_init(&cl_final, (size_t)N, sp.L, sp.rc);
        cell_list_build(&cl_final, p_final, sp.L);

        NeighborList nl_final;
        nbl_init(&nl_final, (size_t)N, sp.rc, 0.3 * sp.rc);
        nbl_build(&nl_final, &cl_final, p_final, sp.L, sp.rc, (size_t)N);

        md_compute_forces_nbl(p_final, &sp, &nl_final);

        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_nbl(p_final, &sp, &nl_final, &cl_final, &K);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }

        free(cl_final.counts);
        free(cl_final.cells);
        free(nl_final.nb);
        free(nl_final.nb_index);
        free(nl_final.prev);
    }
    else if (fastest_method == 3) {
        CellList cl_mt_final;
        cell_list_init(&cl_mt_final, (size_t)N, sp.L, sp.rc);
        cell_list_build(&cl_mt_final, p_final, sp.L);

        md_compute_forces_cell_mt(p_final, &sp, &cl_mt_final);

        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_cell_mt(p_final, &sp, &cl_mt_final, &K);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }

        free(cl_mt_final.counts);
        free(cl_mt_final.cells);
    }
    else if (fastest_method == 4) {
        CellList cl_final;
        cell_list_init(&cl_final, (size_t)N, sp.L, sp.rc);
        cell_list_build(&cl_final, p_final, sp.L);

        NeighborList nl_final;
        nbl_init(&nl_final, (size_t)N, sp.rc, 0.3 * sp.rc);
        nbl_build_pthread(&nl_final, &cl_final, p_final, sp.L, sp.rc, (size_t)N);

        md_compute_forces_nbl_pthread(p_final, &sp, &nl_final);

        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_nbl_pthread(p_final, &sp, &nl_final, &cl_final, &K);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }

        free(cl_final.counts);
        free(cl_final.cells);
        free(nl_final.nb);
        free(nl_final.nb_index);
        free(nl_final.prev);
    }
    else if (fastest_method == 5) {
        CellList cl_final;
        cell_list_init(&cl_final, (size_t)N, sp.L, sp.rc);
        cell_list_build(&cl_final, p_final, sp.L);

        md_compute_forces_cell_manhattan_pthread(p_final, &sp, &cl_final, sp.nthreads);

        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_manhattan_pthread(p_final, &sp, &cl_final, &K, sp.nthreads);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }

        free(cl_final.counts);
        free(cl_final.cells);
    }
    else if (fastest_method == 6) {
        // GPU Manhattan: GPU builds Manhattan cell bins internally each step.
        md_compute_forces_cell_manhattan_gpu(p_final, &sp);

        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_manhattan_gpu(p_final, &sp, &K);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_final = interval(time_start, time_stop);

    fclose(f_final);
    io_write_pdb("output/final_positions.pdb", p_final, (size_t)N);

    printf("Final simulation complete. Time: %.6f seconds\n", time_final);
    printf("Output: output/final_energies.csv, output/final_positions.pdb\n");

    // ------------------------------------------------
    // Cleanup
    // ------------------------------------------------
    nbl_pthread_cleanup();

    free(cl.counts);   free(cl.cells);
    free(cl2.counts);  free(cl2.cells);

    free(nl.nb); free(nl.nb_index); free(nl.prev);

    free(p0_dyn); free(p_full); free(p_cell); free(p_nbl);
    free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);

    return 0;
}