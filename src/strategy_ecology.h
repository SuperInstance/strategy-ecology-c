#ifndef STRATEGY_ECOLOGY_H
#define STRATEGY_ECOLOGY_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Species ──────────────────────────────────────────────────────── */

typedef enum {
    SPECIES_EXPLORER    = 0,
    SPECIES_DIPLOMAT    = 1,
    SPECIES_MARKSMAN    = 2,
    SPECIES_CLIMBER     = 3,
    SPECIES_PROSPECTOR  = 4,
    SPECIES_COUNT       = 5
} Species;

typedef struct {
    double win_rate;   /* typical win rate [0,1] */
    double entropy;    /* typical Shannon entropy of action distribution */
    const char *name;
} SpeciesTraits;

/* Get static traits for a species. */
const SpeciesTraits *species_traits(Species s);

/* ── Population ───────────────────────────────────────────────────── */

typedef struct {
    unsigned int counts[SPECIES_COUNT];
} Population;

Population population_create(const unsigned int counts[SPECIES_COUNT]);
unsigned int population_total(const Population *p);
double      population_fraction(const Population *p, Species s);

/* Shannon diversity index H = -Σ p_i ln(p_i) */
double population_shannon(const Population *p);

/* Simpson diversity index D = 1 - Σ p_i² */
double population_simpson(const Population *p);

/* ── Lotka-Volterra dynamics ──────────────────────────────────────── */

typedef struct {
    double alpha[SPECIES_COUNT][SPECIES_COUNT]; /* interaction matrix */
    double growth[SPECIES_COUNT];               /* intrinsic growth rates */
    double dt;                                  /* time step */
} LotkaVolterra;

/* Build a LotkaVolterra model with default interaction matrix. */
LotkaVolterra lotka_volterra_default(double dt);

/* Run one Euler step. Returns updated population (rounded to ints). */
Population lotka_volterra_step(const LotkaVolterra *lv,
                               const Population *pop);

/* Run n Euler steps. */
Population lotka_volterra_simulate(const LotkaVolterra *lv,
                                   const Population *pop,
                                   size_t steps);

/* ── Ecological Resilience ────────────────────────────────────────── */

typedef struct {
    int surviving[SPECIES_COUNT];   /* 1 if species count > 0 */
    int surviving_count;
    double resilience_index;        /* fraction of species surviving, weighted by evenness */
} EcologicalResilience;

EcologicalResilience ecological_resilience(const Population *before,
                                           const Population *after);

/* ── Species Classification ───────────────────────────────────────── */

typedef struct {
    double win_rate;
    double entropy;
    double cooperation_rate;  /* fraction of cooperative actions */
    double risk_tolerance;    /* [0,1] */
} AgentBehavior;

/* Classify an agent into a species based on behavior profile. */
Species species_classify(const AgentBehavior *behavior);

#ifdef __cplusplus
}
#endif

#endif /* STRATEGY_ECOLOGY_H */
