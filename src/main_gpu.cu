#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "md.h"
#include "gpu_manhattan.h"

static double wrap0L(double x, double L) {
  x = fmod(x, L);
  if (x < 0) x += L;
  return x;
}

static double urand() { return (double)rand() / (double)RAND_MAX; }

static void init_lattice_jitter(Particle *p, size_t N, double L, double jitter_frac) {
  size_t n3 = (size_t)ceil(cbrt((double)N));
  double spacing = L / (double)n3;
  size_t idx = 0;

  for (size_t iz = 0; iz < n3 && idx < N; iz++)
    for (size_t iy = 0; iy < n3 && idx < N; iy++)
      for (size_t ix = 0; ix < n3 && idx < N; ix++) {

        double x = (ix + 0.5) * spacing;
        double y = (iy + 0.5) * spacing;
        double z = (iz + 0.5) * spacing;

        double j = jitter_frac * spacing;
        x += (urand() - 0.5) * 2.0 * j;
        y += (urand() - 0.5) * 2.0 * j;
        z += (urand() - 0.5) * 2.0 * j;

        p[idx].x = wrap0L(x, L);
        p[idx].y = wrap0L(y, L);
        p[idx].z = wrap0L(z, L);

        p[idx].vx = p[idx].vy = p[idx].vz = 0.0;
        p[idx].fx = p[idx].fy = p[idx].fz = 0.0;
        idx++;
      }
}

int main(int argc, char **argv) {
  size_t N = 200000;
  double L = 300.0;

  if (argc >= 2) N = (size_t)atoll(argv[1]);
  if (argc >= 3) L = atof(argv[2]);

  Particle *p = (Particle*)malloc(N * sizeof(Particle));
  if (!p) {
    fprintf(stderr, "malloc failed\n");
    return 1;
  }

  SimParams sp;
  sp.N = N;
  sp.L = L;
  sp.rc = 2.5;
  sp.rc2 = sp.rc * sp.rc;
  sp.sigma = 1.0;
  sp.epsilon = 1.0;
  sp.dt = 0.005;
  sp.nthreads = 1;

  init_lattice_jitter(p, N, L, 0.05);

  double U = md_compute_forces_cell_manhattan_gpu(p, &sp);

  double fsum = 0.0;
  for (size_t i = 0; i < (N < 1000 ? N : 1000); i++) {
    fsum += p[i].fx + p[i].fy + p[i].fz;
  }

  printf("GPU Manhattan test: N=%zu L=%.1f U=%0.6e fsum(first<=1000)=%0.6e\n",
         N, L, U, fsum);

  free(p);
  return 0;
}