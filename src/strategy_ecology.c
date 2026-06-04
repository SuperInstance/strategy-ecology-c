#include "strategy_ecology.h"
#include <math.h>
#include <string.h>
#include <float.h>

/* ── Species traits (static) ──────────────────────────────────── */

static const SpeciesTraits traits_table[SPECIES_COUNT] = {
    /* SPECIES_EXPLORER   */ { 0.45, 2.80, "Explorer"    },
    /* SPECIES_DIPLOMAT   */ { 0.50, 1.50, "Diplomat"    },
    /* SPECIES_MARKSMAN   */ { 0.65, 0.80, "Marksman"    },
    /* SPECIES_CLIMBER    */ { 0.40, 2.20, "Climber"     },
    /* SPECIES_PROSPECTOR */ { 0.55, 1.90, "Prospector"  },
};

const SpeciesTraits *species_traits(Species s) {
    if (s < 0 || s >= SPECIES_COUNT) return NULL;
    return &traits_table[s];
}

/* ── Population ───────────────────────────────────────────────── */

Population population_create(const unsigned int counts[SPECIES_COUNT]) {
    Population p;
    memcpy(p.counts, counts, sizeof(p.counts));
    return p;
}

unsigned int population_total(const Population *p) {
    unsigned int total = 0;
    for (int i = 0; i < SPECIES_COUNT; i++) total += p->counts[i];
    return total;
}

double population_fraction(const Population *p, Species s) {
    unsigned int total = population_total(p);
    if (total == 0 || s < 0 || s >= SPECIES_COUNT) return 0.0;
    return (double)p->counts[s] / total;
}

double population_shannon(const Population *p) {
    unsigned int total = population_total(p);
    if (total == 0) return 0.0;
    double h = 0.0;
    for (int i = 0; i < SPECIES_COUNT; i++) {
        if (p->counts[i] == 0) continue;
        double pi = (double)p->counts[i] / total;
        h -= pi * log(pi);
    }
    return h;
}

double population_simpson(const Population *p) {
    unsigned int total = population_total(p);
    if (total == 0) return 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < SPECIES_COUNT; i++) {
        double pi = (double)p->counts[i] / total;
        sum_sq += pi * pi;
    }
    return 1.0 - sum_sq;
}

/* ── Lotka-Volterra ───────────────────────────────────────────── */

LotkaVolterra lotka_volterra_default(double dt) {
    LotkaVolterra lv;
    lv.dt = dt;

    /* Intrinsic growth rates */
    lv.growth[SPECIES_EXPLORER]   =  0.10;
    lv.growth[SPECIES_DIPLOMAT]   =  0.08;
    lv.growth[SPECIES_MARKSMAN]   =  0.12;
    lv.growth[SPECIES_CLIMBER]    =  0.07;
    lv.growth[SPECIES_PROSPECTOR] =  0.09;

    /* Default interaction matrix: competition model
       alpha[i][j] = effect of species j on species i
       Diagonal = intraspecific competition (self-limiting)
       Off-diag = interspecific competition (mild) */
    static const double default_alpha[SPECIES_COUNT][SPECIES_COUNT] = {
        /* E    D    M    C    P  */
        { 0.02, 0.01, 0.005, 0.008, 0.01 },  /* Explorer */
        { 0.008, 0.02, 0.01, 0.005, 0.012 }, /* Diplomat */
        { 0.005, 0.01, 0.025, 0.01, 0.008 }, /* Marksman */
        { 0.01, 0.008, 0.005, 0.02, 0.01 },  /* Climber */
        { 0.012, 0.01, 0.008, 0.01, 0.02 },  /* Prospector */
    };
    memcpy(lv.alpha, default_alpha, sizeof(default_alpha));
    return lv;
}

Population lotka_volterra_step(const LotkaVolterra *lv,
                               const Population *pop) {
    double n[SPECIES_COUNT];
    for (int i = 0; i < SPECIES_COUNT; i++) {
        n[i] = (double)pop->counts[i];
    }

    /* Euler step: dn_i/dt = r_i * n_i * (1 - Σ_j α_ij * n_j) */
    double dn[SPECIES_COUNT];
    for (int i = 0; i < SPECIES_COUNT; i++) {
        double interaction = 0.0;
        for (int j = 0; j < SPECIES_COUNT; j++) {
            interaction += lv->alpha[i][j] * n[j];
        }
        dn[i] = lv->growth[i] * n[i] * (1.0 - interaction);
    }

    Population next;
    for (int i = 0; i < SPECIES_COUNT; i++) {
        double val = n[i] + dn[i] * lv->dt;
        if (val < 0.0) val = 0.0;
        next.counts[i] = (unsigned int)(val + 0.5); /* round */
    }
    return next;
}

Population lotka_volterra_simulate(const LotkaVolterra *lv,
                                   const Population *pop,
                                   size_t steps) {
    Population current = *pop;
    for (size_t s = 0; s < steps; s++) {
        current = lotka_volterra_step(lv, &current);
    }
    return current;
}

/* ── Ecological Resilience ────────────────────────────────────── */

EcologicalResilience ecological_resilience(const Population *before,
                                           const Population *after) {
    EcologicalResilience er;
    er.surviving_count = 0;

    for (int i = 0; i < SPECIES_COUNT; i++) {
        er.surviving[i] = (after->counts[i] > 0) ? 1 : 0;
        er.surviving_count += er.surviving[i];
    }

    /* Resilience index: fraction surviving * evenness factor */
    double frac_surviving = (double)er.surviving_count / SPECIES_COUNT;

    /* Use evenness from after-population (0 if empty) */
    double evenness = 0.0;
    unsigned int total = population_total(after);
    if (total > 0) {
        double h = population_shannon(after);
        evenness = h / log((double)SPECIES_COUNT);
    }

    er.resilience_index = frac_surviving * (0.7 + 0.3 * evenness);
    return er;
}

/* ── Species Classification ───────────────────────────────────── */

Species species_classify(const AgentBehavior *b) {
    /* Score each species by matching behavior to typical traits.
       We combine win_rate, entropy, cooperation, and risk tolerance. */

    /* Expected profiles per species:
       Explorer:   low-mid win, high entropy, mid coop, high risk
       Diplomat:   mid win, low entropy, high coop, low risk
       Marksman:   high win, low entropy, low coop, mid risk
       Climber:    low win, mid-high entropy, low coop, high risk
       Prospector: mid-high win, mid entropy, mid coop, mid-high risk */

    double scores[SPECIES_COUNT];

    /* Explorer */
    scores[SPECIES_EXPLORER] =
        (1.0 - fabs(b->win_rate - 0.45)) +
        (1.0 - fabs(b->entropy - 2.80) / 3.0) +
        (1.0 - fabs(b->cooperation_rate - 0.4)) +
        (b->risk_tolerance * 1.5);

    /* Diplomat */
    scores[SPECIES_DIPLOMAT] =
        (1.0 - fabs(b->win_rate - 0.50)) +
        (1.0 - fabs(b->entropy - 1.50) / 3.0) +
        (b->cooperation_rate * 1.5) +
        ((1.0 - b->risk_tolerance) * 1.2);

    /* Marksman */
    scores[SPECIES_MARKSMAN] =
        (b->win_rate * 1.5) +
        (1.0 - fabs(b->entropy - 0.80) / 3.0) +
        ((1.0 - b->cooperation_rate) * 1.2) +
        (1.0 - fabs(b->risk_tolerance - 0.5));

    /* Climber */
    scores[SPECIES_CLIMBER] =
        (1.0 - fabs(b->win_rate - 0.40)) +
        (1.0 - fabs(b->entropy - 2.20) / 3.0) +
        ((1.0 - b->cooperation_rate) * 1.0) +
        (b->risk_tolerance * 1.8);

    /* Prospector */
    scores[SPECIES_PROSPECTOR] =
        (1.0 - fabs(b->win_rate - 0.55)) +
        (1.0 - fabs(b->entropy - 1.90) / 3.0) +
        (1.0 - fabs(b->cooperation_rate - 0.5)) +
        (b->risk_tolerance * 1.3);

    /* Pick highest score */
    Species best = SPECIES_EXPLORER;
    for (int i = 1; i < SPECIES_COUNT; i++) {
        if (scores[i] > scores[best]) best = (Species)i;
    }
    return best;
}
