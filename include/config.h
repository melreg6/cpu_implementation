#ifndef CONFIG_H
#define CONFIG_H

// -----------------------------
// Lennard-Jones (unitless)
// Positions are in Angstroms but σ and ε define dimensionless LJ units
// -----------------------------
#define CONF_SIGMA       1.0     // LJ sigma (defines distance scale)
#define CONF_EPSILON     1.0     // LJ epsilon (defines energy scale)

// -----------------------------
// LJ Cutoff (in Angstrom units)
// -----------------------------
#define CONF_CUTOFF      30

// -----------------------------
// Time step + steps
// -----------------------------
#define CONF_DELTA_T     0.001
#define AUTOTUNE_N_TIMESTEPS 2000

#define USER_N_TIMESTEPS 4000

// -----------------------------
// Box min/max (fallback — will be overwritten by PDB bounds)
// -----------------------------
#define CONF_BOX_MIN     0.0
#define CONF_BOX_MAX     100.0

#define N_THREADS        16

#endif
