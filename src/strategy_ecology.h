#ifndef STRATEGY_ECOLOGY_H
#define STRATEGY_ECOLOGY_H

#include <stddef.h>
#include <stdbool.h>

/*
 * Strategy Ecology — Lotka-Volterra style species dynamics for ternary agents.
 *
 * Each agent carries a ternary strategy in {0, 1, 2}.
 * Agents are classified into 5 species categories based on strategy and fitness.
 * Species interact via a Lotka-Volterra system with Euler integration.
 */

/* Ternary strategy values */
typedef enum { STRAT_ROCK = 0, STRAT_PAPER = 1, STRAT_SCISSORS = 2 } TernaryStrategy;

/* 5 species categories */
typedef enum {
    SPECIES_ALPHA   = 0,  /* dominant strategy */
    SPECIES_BETA    = 1,  /* counter to dominant */
    SPECIES_GAMMA   = 2,  /* counter to counter */
    SPECIES_DELTA   = 3,  /* neutral / generalist */
    SPECIES_EPSILON = 4,  /* low-fitness scavenger */
    SPECIES_COUNT   = 5
} SpeciesID;

/* Agent */
typedef struct {
    TernaryStrategy strategy;
    double fitness;
    SpeciesID species;
} Agent;

/* Population of N agents */
typedef struct {
    Agent *agents;
    size_t count;
    double payoff_matrix[3][3];  /* payoff[a][b] = payoff to strategy a facing b */
} Population;

/* Species tracker — counts per species */
typedef struct {
    int counts[SPECIES_COUNT];
    double total_fitness[SPECIES_COUNT];
} SpeciesTracker;

/* Lotka-Volterra parameters and state */
typedef struct {
    double populations[SPECIES_COUNT];   /* current population sizes */
    double growth_rates[SPECIES_COUNT];  /* intrinsic growth rates r_i */
    double interaction[SPECIES_COUNT][SPECIES_COUNT]; /* A_ij interaction matrix */
    double carrying_capacity[SPECIES_COUNT];          /* K_i */
    double dt;                                          /* time step */
    double decay_rate;                                  /* universal decay */
    int generations;                                    /* number of integration steps */
} LotkaVolterra;

/* Coexistence check result */
typedef struct {
    bool all_survive;
    double min_population;
    int surviving_count;
    double final_populations[SPECIES_COUNT];
} CoexistenceResult;

/* Fitness landscape — maps strategy + context to fitness */
typedef struct {
    double peaks[3];         /* fitness peak per strategy */
    double widths[3];        /* width/breadth of each peak */
    double shift;            /* environmental shift factor */
    double noise;            /* noise amplitude */
    int seed;                /* RNG seed for reproducibility */
} FitnessLandscape;

/* --- Population API --- */

/* Create a population with N random agents and default payoff matrix */
Population *population_create(size_t n, int seed);

/* Free a population */
void population_free(Population *pop);

/* Evaluate all pairwise interactions, updating agent fitness */
void population_interact(Population *pop);

/* Classify each agent into a species based on strategy and fitness */
void population_classify(Population *pop);

/* --- SpeciesTracker API --- */

/* Initialize a tracker from a classified population */
void species_tracker_init(SpeciesTracker *st, const Population *pop);

/* Get count for a species */
int species_tracker_count(const SpeciesTracker *st, SpeciesID id);

/* Get average fitness for a species; returns 0 if count is 0 */
double species_tracker_avg_fitness(const SpeciesTracker *st, SpeciesID id);

/* --- LotkaVolterra API --- */

/* Initialize with default competitive Lotka-Volterra parameters */
void lotka_volterra_init(LotkaVolterra *lv, const SpeciesTracker *st);

/* Perform one Euler integration step */
void lotka_volterra_step(LotkaVolterra *lv);

/* Run for N generations, updating populations in place */
void lotka_volterra_run(LotkaVolterra *lv, int generations);

/* --- CoexistenceCheck API --- */

/* Check if all species survive after N generations */
CoexistenceResult coexistence_check(const LotkaVolterra *lv, int generations, double threshold);

/* --- FitnessLandscape API --- */

/* Initialize a fitness landscape */
void fitness_landscape_init(FitnessLandscape *fl, double shift, double noise, int seed);

/* Evaluate fitness for a given strategy */
double fitness_landscape_evaluate(const FitnessLandscape *fl, TernaryStrategy s);

/* Apply fitness landscape to a population */
void fitness_landscape_apply(const FitnessLandscape *fl, Population *pop);

#endif /* STRATEGY_ECOLOGY_H */
