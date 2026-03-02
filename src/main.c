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

#define MAX_ATOMS 50000

// ----------------------------------------------------
// Time measurement by clock_gettime()
// ----------------------------------------------------
double interval(struct timespec start, struct timespec end)
{
    struct timespec temp;
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    if (temp.tv_nsec < 0) {
        temp.tv_sec = temp.tv_sec - 1;
        temp.tv_nsec = temp.tv_nsec + 1000000000;
    }
    return (((double)temp.tv_sec) + ((double)temp.tv_nsec) * 1.0e-9);
}

// ----------------------------------------------------
// Helper: copy particle array
// ----------------------------------------------------
static void copy_particles(Particle *dst, const Particle *src, size_t N)
{
    for (size_t i = 0; i < N; i++) {
        dst[i] = src[i];
    }
}

// ----------------------------------------------------
// Helper: open energy CSV file
// ----------------------------------------------------
static FILE *open_energy_csv(const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) {
        perror("fopen energy csv");
        return NULL;
    }
    fprintf(f, "step,K,U,E\n");
    return f;
}

// ----------------------------------------------------
// Helper: wrap position into [0, L)
// ----------------------------------------------------
static inline double wrap_pos(double x, double L)
{
    x = fmod(x, L);
    if (x < 0) x += L;
    return x;
}

// ----------------------------------------------------
// Helper: kinetic energy (assumes unit mass)
// ----------------------------------------------------
static double kinetic_energy(const Particle *p, size_t N)
{
    double K = 0.0;
    for (size_t i = 0; i < N; i++) {
        K += 0.5 * (p[i].vx * p[i].vx + p[i].vy * p[i].vy + p[i].vz * p[i].vz);
    }
    return K;
}

// ----------------------------------------------------
// Manhattan Velocity-Verlet integrator (Pthread force)
// ----------------------------------------------------
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

// ----------------------------------------------------
// main
// ----------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage: %s input_path [--threads N | -t N]\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];

    printf("\n================ INITIALIZING SIMULATION ================\n");

    // ---- parse optional --threads / -t ----
    int user_threads = N_THREADS;
    for (int i = 2; i < argc; i++) {
        if ((!strcmp(argv[i], "--threads") || !strcmp(argv[i], "-t")) && (i + 1 < argc)) {
            user_threads = atoi(argv[++i]);
            if (user_threads < 1) user_threads = 1;
        }
    }

    double xmin, xmax, ymin, ymax, zmin, zmax;

    printf("Converting input file via Python script...\n");
    char command[256];
    sprintf(command, "./venv/bin/python3 src/py/converter.py %s temp_coords.bin temp_metadata.txt", input_path);
    int status = system(command);

    if (status != 0) {
        fprintf(stderr, "Error: Python conversion failed.\n");
        return 1;
    }

    FILE *f_meta = fopen("temp_metadata.txt", "r");
    int N;
    if (f_meta) {
        fscanf(f_meta, "%d", &N);
        fclose(f_meta);
    } else {
        fprintf(stderr, "Metadata not found!\n");
        return 1;
    }

    Particle *p0 = malloc((size_t)N * sizeof(Particle));
    if (!p0) {
        fprintf(stderr, "Allocation failure p0\n");
        return 1;
    }

    int foo;
    foo = binary_importer("temp_coords.bin", p0, N,
                         &xmin, &xmax,
                         &ymin, &ymax,
                         &zmin, &zmax);
    (void)foo;

    if (N <= 0) {
        fprintf(stderr, "Failed to read PDB: %s\n", input_path);
        free(p0);
        return 1;
    }

    printf("deleting temp files \n");
    (void)system("rm -rf temp_coords.bin");
    (void)system("rm -rf temp_metadata.txt");

    printf("Loaded %d atoms from %s\n", N, input_path);
    printf("Bounds: dx = %.3f  dy = %.3f  dz = %.3f\n",
           xmax - xmin, ymax - ymin, zmax - zmin);

    // ------------------------------------------------
    // 2. Set up simulation parameters
    // ------------------------------------------------
    SimParams sp;
    sp.N        = (size_t)N;
    sp.sigma    = CONF_SIGMA;
    sp.epsilon  = CONF_EPSILON;
    sp.rc       = CONF_CUTOFF;
    sp.rc2      = sp.rc * sp.rc;
    sp.dt       = CONF_DELTA_T;

    // Use CLI override
    sp.nthreads = user_threads;

    // Use largest span as box length
    sp.L = fmax(fmax(xmax - xmin, ymax - ymin), zmax - zmin);
    if (sp.L <= 0.0) sp.L = CONF_BOX_MAX;

    printf("Simulation box L = %.3f\n", sp.L);
    printf("rc = %.3f, dt = %.5f, steps = %d\n",
           sp.rc, sp.dt, AUTOTUNE_N_TIMESTEPS);
    printf("Manual thread override: %d\n", sp.nthreads);

    // ------------------------------------------------
    // 3. Make copies for each MD method and for the final simulation
    // ------------------------------------------------
    Particle *p_full = malloc((size_t)N * sizeof(Particle));
    Particle *p_cell = malloc((size_t)N * sizeof(Particle));
    Particle *p_nbl  = malloc((size_t)N * sizeof(Particle));
    Particle *p_cell_mt = malloc((size_t)N * sizeof(Particle));
    Particle *p_nbl_pthread = malloc((size_t)N * sizeof(Particle));
    Particle *p_manhattan = malloc((size_t)N * sizeof(Particle));
    Particle *p_final = malloc((size_t)N * sizeof(Particle));

    if (!p_full || !p_cell || !p_nbl || !p_cell_mt || !p_nbl_pthread || !p_manhattan || !p_final) {
        fprintf(stderr, "Allocation failure for particle copies\n");
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    copy_particles(p_full, p0, (size_t)N);
    copy_particles(p_cell, p0, (size_t)N);
    copy_particles(p_nbl,  p0, (size_t)N);
    copy_particles(p_cell_mt, p0, (size_t)N);
    copy_particles(p_nbl_pthread, p0, (size_t)N);
    copy_particles(p_manhattan, p0, (size_t)N);
    copy_particles(p_final, p0, (size_t)N);

    printf("\n================ BEGINNING AUTOTUNING ================\n");

    struct timespec time_start, time_stop;

    // ------------------------------------------------
    // 4. FULL O(N^2) MD - WITH TIMING
    // ------------------------------------------------
    printf("\n================ FULL O(N^2) MD ================\n");

    FILE *f_full = open_energy_csv("output/full_energies.csv");
    if (!f_full) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    (void) md_compute_forces_full(p_full, &sp);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_full(p_full, &sp, &K);
        fprintf(f_full, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_full = interval(time_start, time_stop);

    fclose(f_full);
    io_write_pdb("output/full_positions.pdb", p_full, (size_t)N);
    printf("FULL MD done. Time: %.6f seconds\n", time_full);

    // ------------------------------------------------
    // 5. CELL-LIST MD - WITH TIMING
    // ------------------------------------------------
    printf("\n================ CELL-LIST MD ================\n");

    FILE *f_cell = open_energy_csv("output/cell_energies.csv");
    if (!f_cell) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    CellList cl;
    cell_list_init(&cl, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl, p_cell, sp.L);
    (void) md_compute_forces_cell(p_cell, &sp, &cl);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_cell(p_cell, &sp, &cl, &K);
        fprintf(f_cell, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_cell = interval(time_start, time_stop);

    fclose(f_cell);
    io_write_pdb("output/cell_positions.pdb", p_cell, (size_t)N);
    printf("CELL-LIST MD done. Time: %.6f seconds\n", time_cell);

    // ------------------------------------------------
    // 6. NEIGHBOR-LIST MD (Serial) - WITH TIMING
    // ------------------------------------------------
    printf("\n================ NEIGHBOR-LIST MD (Serial) ================\n");

    FILE *f_nbl = open_energy_csv("output/nbl_energies.csv");
    if (!f_nbl) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    CellList cl2;
    cell_list_init(&cl2, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl2, p_nbl, sp.L);

    NeighborList nl;
    nbl_init(&nl, (size_t)N, sp.rc, 0.3 * sp.rc);
    nbl_build(&nl, &cl2, p_nbl, sp.L, sp.rc, (size_t)N);

    (void) md_compute_forces_nbl(p_nbl, &sp, &nl);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_nbl(p_nbl, &sp, &nl, &cl2, &K);
        fprintf(f_nbl, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_nbl = interval(time_start, time_stop);

    fclose(f_nbl);
    io_write_pdb("output/nbl_positions.pdb", p_nbl, (size_t)N);
    printf("NEIGHBOR-LIST MD (Serial) done. Time: %.6f seconds\n", time_nbl);

    // ------------------------------------------------
    // HALF-SHELL MULTI-THREADED CELL-LIST MD
    // ------------------------------------------------
    printf("\n================ HALF-SHELL PARALLEL CELL-LIST MD ================\n");

    FILE *f_cell_mt = open_energy_csv("output/cell_half_mt_energies.csv");
    if (!f_cell_mt){
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    CellList cl_half;
    cell_list_init(&cl_half, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl_half, p_cell_mt, sp.L);

    md_compute_forces_cell_mt(p_cell_mt, &sp, &cl_half);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_cell_mt(p_cell_mt, &sp, &cl_half, &K);
        fprintf(f_cell_mt, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_cell_mt = interval(time_start, time_stop);

    fclose(f_cell_mt);
    io_write_pdb("output/cell_half_mt_positions.pdb", p_cell_mt, (size_t)N);
    printf("HALF-SHELL MULTI-THREADED CELL-LIST MD done. Time: %.6f seconds\n", time_cell_mt);

    // ------------------------------------------------
    // 7. NEIGHBOR-LIST MD (Pthread) - WITH TIMING
    // ------------------------------------------------
    printf("\n================ NEIGHBOR-LIST MD (Pthread) ================\n");
    nbl_pthread_print_info();

    FILE *f_nbl_pthread = open_energy_csv("output/nbl_pthread_energies.csv");
    if (!f_nbl_pthread) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        free(p_cell_mt); free(p_nbl_pthread); free(p_manhattan); free(p_final);
        return 1;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    nbl_pthread_init((size_t)N);

    CellList cl3;
    cell_list_init(&cl3, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl3, p_nbl_pthread, sp.L);

    NeighborList nl_pthread;
    nbl_init(&nl_pthread, (size_t)N, sp.rc, 0.3 * sp.rc);
    nbl_build_pthread(&nl_pthread, &cl3, p_nbl_pthread, sp.L, sp.rc, (size_t)N);

    (void) md_compute_forces_nbl_pthread(p_nbl_pthread, &sp, &nl_pthread);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_nbl_pthread(p_nbl_pthread, &sp, &nl_pthread, &cl3, &K);
        fprintf(f_nbl_pthread, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_nbl_pthread = interval(time_start, time_stop);

    fclose(f_nbl_pthread);
    io_write_pdb("output/nbl_pthread_positions.pdb", p_nbl_pthread, (size_t)N);
    printf("NEIGHBOR-LIST MD (Pthread) done. Time: %.6f seconds\n", time_nbl_pthread);
    printf("Memory usage for private arrays: %.2f MB\n",
           nbl_pthread_memory_usage((size_t)N) / (1024.0 * 1024.0));

    // ------------------------------------------------
    // 7.5 MANHATTAN CELL-LIST MD (Pthread) - WITH TIMING
    // ------------------------------------------------
    printf("\n================ MANHATTAN CELL-LIST MD (Pthread) ================\n");
    printf("Manhattan threads: %d\n", sp.nthreads);

    FILE *f_manh = open_energy_csv("output/manhattan_pthread_energies.csv");
    if (!f_manh) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
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

    printf("\nFASTEST METHOD: %s (%.6f seconds)\n", fastest_name, fastest_time);

    // ------------------------------------------------
    // 9. RUN FULL SIMULATION WITH FASTEST METHOD
    // ------------------------------------------------
    printf("\n========================================\n");
    printf("RUNNING FULL SIMULATION WITH: %s\nfor %d TIMESTEPS\n",
           fastest_name, USER_N_TIMESTEPS);
    printf("========================================\n");

    FILE *f_final = open_energy_csv("output/final_energies.csv");
    if (!f_final) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
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
    free(cl3.counts);  free(cl3.cells);
    free(cl_half.counts); free(cl_half.cells);
    free(cl_manh.counts); free(cl_manh.cells);

    free(nl.nb); free(nl.nb_index); free(nl.prev);
    free(nl_pthread.nb); free(nl_pthread.nb_index); free(nl_pthread.prev);

    free(p0);
    free(p_full);
    free(p_cell);
    free(p_nbl);
    free(p_cell_mt);
    free(p_nbl_pthread);
    free(p_manhattan);
    free(p_final);

    printf("\nAll simulations complete.\n");
    return 0;
}