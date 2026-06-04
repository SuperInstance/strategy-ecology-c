#include "strategy_ecology.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Simple LCG RNG for reproducibility */
static unsigned int rng_state;

static unsigned int rng_next(void) {
    rng_state = rng_state * 1103515245u + 12345u;
    return rng_state;
}

static double rng_uniform(void) {
    return (double)(rng_next() & 0x7fffffff) / (double)0x7fffffff;
}

/* --- Population API --- */

Population *population_create(size_t n, int seed) {
    Population *pop = (Population *)calloc(1, sizeof(Population));
    if (!pop) return NULL;

    pop->count = n;
    pop->agents = (Agent *)calloc(n, sizeof(Agent));
    if (!pop->agents) { free(pop); return NULL; }

    rng_state = (unsigned int)seed;

    /* Default RPS payoff matrix:
     * Rock vs Rock=0.5, Rock vs Paper=0, Rock vs Scissors=1
     * Paper vs Rock=1, Paper vs Paper=0.5, Paper vs Scissors=0
     * Scissors vs Rock=0, Scissors vs Paper=1, Scissors vs Scissors=0.5
     */
    double pm[3][3] = {
        {0.5, 0.0, 1.0},
        {1.0, 0.5, 0.0},
        {0.0, 1.0, 0.5}
    };
    memcpy(pop->payoff_matrix, pm, sizeof(pm));

    for (size_t i = 0; i < n; i++) {
        pop->agents[i].strategy = (TernaryStrategy)(rng_next() % 3);
        pop->agents[i].fitness = 0.0;
        pop->agents[i].species = SPECIES_ALPHA;
    }

    return pop;
}

void population_free(Population *pop) {
    if (pop) {
        free(pop->agents);
        free(pop);
    }
}

void population_interact(Population *pop) {
    if (!pop || pop->count == 0) return;

    /* Reset fitness */
    for (size_t i = 0; i < pop->count; i++) {
        pop->agents[i].fitness = 0.0;
    }

    /* All-pairs interaction */
    for (size_t i = 0; i < pop->count; i++) {
        for (size_t j = i + 1; j < pop->count; j++) {
            int si = pop->agents[i].strategy;
            int sj = pop->agents[j].strategy;
            double payoff_i = pop->payoff_matrix[si][sj];
            double payoff_j = pop->payoff_matrix[sj][si];
            pop->agents[i].fitness += payoff_i;
            pop->agents[j].fitness += payoff_j;
        }
    }

    /* Normalize by number of interactions */
    if (pop->count > 1) {
        double norm = 1.0 / (double)(pop->count - 1);
        for (size_t i = 0; i < pop->count; i++) {
            pop->agents[i].fitness *= norm;
        }
    }
}

void population_classify(Population *pop) {
    if (!pop || pop->count == 0) return;

    /* First pass: find fitness range */
    double min_fit = pop->agents[0].fitness;
    double max_fit = pop->agents[0].fitness;
    for (size_t i = 1; i < pop->count; i++) {
        if (pop->agents[i].fitness < min_fit) min_fit = pop->agents[i].fitness;
        if (pop->agents[i].fitness > max_fit) max_fit = pop->agents[i].fitness;
    }
    double range = max_fit - min_fit;
    double threshold = (range > 1e-12) ? range / 3.0 : 1e-12;

    for (size_t i = 0; i < pop->count; i++) {
        TernaryStrategy s = pop->agents[i].strategy;
        double f = pop->agents[i].fitness;
        double rel = f - min_fit;

        if (rel >= 2.0 * threshold) {
            /* High fitness */
            pop->agents[i].species = (SpeciesID)s;  /* ALPHA, BETA, GAMMA */
        } else if (rel >= threshold) {
            pop->agents[i].species = SPECIES_DELTA;  /* medium = generalist */
        } else {
            pop->agents[i].species = SPECIES_EPSILON; /* low = scavenger */
        }
    }
}

/* --- SpeciesTracker API --- */

void species_tracker_init(SpeciesTracker *st, const Population *pop) {
    memset(st, 0, sizeof(SpeciesTracker));
    if (!pop) return;

    for (size_t i = 0; i < pop->count; i++) {
        SpeciesID sp = pop->agents[i].species;
        st->counts[sp]++;
        st->total_fitness[sp] += pop->agents[i].fitness;
    }
}

int species_tracker_count(const SpeciesTracker *st, SpeciesID id) {
    if (!st || id < 0 || id >= SPECIES_COUNT) return 0;
    return st->counts[id];
}

double species_tracker_avg_fitness(const SpeciesTracker *st, SpeciesID id) {
    if (!st || id < 0 || id >= SPECIES_COUNT || st->counts[id] == 0) return 0.0;
    return st->total_fitness[id] / (double)st->counts[id];
}

/* --- LotkaVolterra API --- */

void lotka_volterra_init(LotkaVolterra *lv, const SpeciesTracker *st) {
    memset(lv, 0, sizeof(LotkaVolterra));

    /* Default parameters: competitive LV with cyclic dominance */
    double r_default[SPECIES_COUNT] = {1.0, 0.9, 0.95, 0.8, 0.7};
    double k_default[SPECIES_COUNT] = {100.0, 100.0, 100.0, 80.0, 60.0};

    /* Interaction matrix: cyclic dominance + competition */
    /* A[i][j] > 0 means j benefits i; < 0 means j harms i */
    double a_default[SPECIES_COUNT][SPECIES_COUNT] = {
        /*  ALP   BET   GAM   DEL   EPS */
        {-0.01,  0.005, -0.008, -0.005,  0.002},  /* ALP: helped by BET, hurt by GAM */
        {-0.008, -0.01,  0.005, -0.005,  0.002},  /* BET: helped by GAM, hurt by ALP */
        { 0.005, -0.008, -0.01,  -0.005,  0.002},  /* GAM: helped by ALP, hurt by BET */
        {-0.003, -0.003, -0.003, -0.01,   0.001},  /* DEL: generalist, mild competition */
        { 0.003,  0.003,  0.003,  0.001, -0.01}    /* EPS: benefits from others' activity */
    };

    memcpy(lv->growth_rates, r_default, sizeof(r_default));
    memcpy(lv->carrying_capacity, k_default, sizeof(k_default));
    memcpy(lv->interaction, a_default, sizeof(a_default));

    lv->dt = 0.01;
    lv->decay_rate = 0.001;
    lv->generations = 0;

    /* Initialize populations from tracker */
    if (st) {
        for (int i = 0; i < SPECIES_COUNT; i++) {
            lv->populations[i] = (double)st->counts[i];
            if (lv->populations[i] < 1.0) lv->populations[i] = 10.0;  /* floor */
        }
    } else {
        for (int i = 0; i < SPECIES_COUNT; i++) {
            lv->populations[i] = 20.0;
        }
    }
}

void lotka_volterra_step(LotkaVolterra *lv) {
    double dxdt[SPECIES_COUNT];

    for (int i = 0; i < SPECIES_COUNT; i++) {
        double N_i = lv->populations[i];
        double r_i = lv->growth_rates[i];
        double K_i = lv->carrying_capacity[i];

        /* Logistic growth term */
        double growth = r_i * N_i * (1.0 - N_i / K_i);

        /* Interaction term: sum_j A_ij * N_i * N_j */
        double interaction_sum = 0.0;
        for (int j = 0; j < SPECIES_COUNT; j++) {
            if (j != i) {
                interaction_sum += lv->interaction[i][j] * N_i * lv->populations[j];
            }
        }

        /* Decay term */
        double decay = lv->decay_rate * N_i;

        dxdt[i] = growth + interaction_sum - decay;
    }

    /* Euler update */
    for (int i = 0; i < SPECIES_COUNT; i++) {
        lv->populations[i] += lv->dt * dxdt[i];
        /* Clamp to non-negative */
        if (lv->populations[i] < 0.0) lv->populations[i] = 0.0;
    }
}

void lotka_volterra_run(LotkaVolterra *lv, int generations) {
    for (int g = 0; g < generations; g++) {
        lotka_volterra_step(lv);
    }
    lv->generations += generations;
}

/* --- CoexistenceCheck API --- */

CoexistenceResult coexistence_check(const LotkaVolterra *lv, int generations, double threshold) {
    CoexistenceResult result;
    result.all_survive = true;
    result.min_population = 1e30;
    result.surviving_count = 0;

    /* Copy state to run simulation */
    LotkaVolterra sim;
    memcpy(&sim, lv, sizeof(LotkaVolterra));

    lotka_volterra_run(&sim, generations);

    for (int i = 0; i < SPECIES_COUNT; i++) {
        result.final_populations[i] = sim.populations[i];
        if (sim.populations[i] >= threshold) {
            result.surviving_count++;
        } else {
            result.all_survive = false;
        }
        if (sim.populations[i] < result.min_population) {
            result.min_population = sim.populations[i];
        }
    }

    return result;
}

/* --- FitnessLandscape API --- */

void fitness_landscape_init(FitnessLandscape *fl, double shift, double noise, int seed) {
    fl->peaks[STRAT_ROCK]     = 1.0;
    fl->peaks[STRAT_PAPER]    = 1.0;
    fl->peaks[STRAT_SCISSORS] = 1.0;

    fl->widths[STRAT_ROCK]     = 0.5;
    fl->widths[STRAT_PAPER]    = 0.5;
    fl->widths[STRAT_SCISSORS] = 0.5;

    fl->shift = shift;
    fl->noise = noise;
    fl->seed  = seed;
    rng_state = (unsigned int)seed;
}

double fitness_landscape_evaluate(const FitnessLandscape *fl, TernaryStrategy s) {
    double base = fl->peaks[s];
    /* Apply shift: strategies respond differently to environmental change */
    double shifted = base - fl->shift * ((double)s - 1.0);
    double gaussian = exp(-(shifted * shifted) / (2.0 * fl->widths[s] * fl->widths[s]));
    /* Add noise */
    double n = fl->noise > 0 ? fl->noise * (2.0 * rng_uniform() - 1.0) : 0.0;
    double fitness = shifted * gaussian + n;
    return fitness > 0 ? fitness : 0.0;
}

void fitness_landscape_apply(const FitnessLandscape *fl, Population *pop) {
    if (!fl || !pop) return;
    rng_state = (unsigned int)fl->seed;
    for (size_t i = 0; i < pop->count; i++) {
        pop->agents[i].fitness = fitness_landscape_evaluate(fl, pop->agents[i].strategy);
    }
}
