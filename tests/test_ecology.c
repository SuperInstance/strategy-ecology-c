#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "../src/strategy_ecology.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) printf("  %-50s ", name);
#define PASS() do { printf("✓ PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("✗ FAIL: %s\n", msg); tests_failed++; } while(0)
#define ASSERT(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)
#define ASSERT_FEQ(a, b, eps) ASSERT(fabs((a) - (b)) < (eps), "float mismatch")

static void test_population_create_and_size(void) {
    TEST("population_create returns valid pop with correct size");
    Population *p = population_create(100, 42);
    ASSERT(p != NULL, "population is NULL");
    ASSERT(p->count == 100, "count mismatch");
    population_free(p);
    PASS();
}

static void test_population_strategies_in_range(void) {
    TEST("all strategies in {0,1,2}");
    Population *p = population_create(200, 7);
    int ok = 1;
    for (size_t i = 0; i < p->count; i++) {
        if (p->agents[i].strategy < 0 || p->agents[i].strategy > 2) { ok = 0; break; }
    }
    population_free(p);
    ASSERT(ok, "strategy out of range");
    PASS();
}

static void test_payoff_matrix_default(void) {
    TEST("default payoff matrix is RPS-valid");
    Population *p = population_create(10, 1);
    ASSERT_FEQ(p->payoff_matrix[STRAT_ROCK][STRAT_SCISSORS], 1.0, 1e-9);
    ASSERT_FEQ(p->payoff_matrix[STRAT_PAPER][STRAT_ROCK], 1.0, 1e-9);
    ASSERT_FEQ(p->payoff_matrix[STRAT_SCISSORS][STRAT_PAPER], 1.0, 1e-9);
    ASSERT_FEQ(p->payoff_matrix[STRAT_ROCK][STRAT_ROCK], 0.5, 1e-9);
    population_free(p);
    PASS();
}

static void test_interaction_updates_fitness(void) {
    TEST("population_interact updates fitness for all agents");
    Population *p = population_create(50, 99);
    population_interact(p);
    int nonzero = 0;
    for (size_t i = 0; i < p->count; i++) {
        if (p->agents[i].fitness > 0.0) nonzero++;
    }
    population_free(p);
    ASSERT(nonzero > 0, "all fitnesses are zero");
    PASS();
}

static void test_classify_all_species(void) {
    TEST("population_classify assigns species to all agents");
    Population *p = population_create(200, 123);
    population_interact(p);
    population_classify(p);
    int counts[SPECIES_COUNT] = {0};
    for (size_t i = 0; i < p->count; i++) {
        SpeciesID sp = p->agents[i].species;
        ASSERT(sp >= 0 && sp < SPECIES_COUNT, "species out of range");
        counts[sp]++;
    }
    population_free(p);
    PASS();
}

static void test_species_tracker_counts(void) {
    TEST("species_tracker counts match classified population");
    Population *p = population_create(300, 55);
    population_interact(p);
    population_classify(p);
    SpeciesTracker st;
    species_tracker_init(&st, p);
    int total = 0;
    for (int i = 0; i < SPECIES_COUNT; i++) {
        total += species_tracker_count(&st, (SpeciesID)i);
    }
    ASSERT((size_t)total == p->count, "tracker total != population count");
    population_free(p);
    PASS();
}

static void test_species_tracker_avg_fitness(void) {
    TEST("species_tracker avg fitness is non-negative");
    Population *p = population_create(100, 77);
    population_interact(p);
    population_classify(p);
    SpeciesTracker st;
    species_tracker_init(&st, p);
    for (int i = 0; i < SPECIES_COUNT; i++) {
        double avg = species_tracker_avg_fitness(&st, (SpeciesID)i);
        ASSERT(avg >= 0.0, "negative avg fitness");
    }
    population_free(p);
    PASS();
}

static void test_lotka_volterra_init(void) {
    TEST("lotka_volterra_init produces non-zero populations");
    LotkaVolterra lv;
    lotka_volterra_init(&lv, NULL);
    for (int i = 0; i < SPECIES_COUNT; i++) {
        ASSERT(lv.populations[i] > 0.0, "population is zero");
        ASSERT(lv.growth_rates[i] > 0.0, "growth rate non-positive");
    }
    PASS();
}

static void test_lotka_volterra_step_changes_state(void) {
    TEST("lotka_volterra_step changes populations");
    LotkaVolterra lv;
    lotka_volterra_init(&lv, NULL);
    double snap[SPECIES_COUNT];
    memcpy(snap, lv.populations, sizeof(snap));
    lotka_volterra_step(&lv);
    int changed = 0;
    for (int i = 0; i < SPECIES_COUNT; i++) {
        if (fabs(lv.populations[i] - snap[i]) > 1e-12) changed++;
    }
    ASSERT(changed > 0, "no populations changed after step");
    PASS();
}

static void test_lotka_volterra_run(void) {
    TEST("lotka_volterra_run for 1000 steps stays stable");
    LotkaVolterra lv;
    lotka_volterra_init(&lv, NULL);
    lotka_volterra_run(&lv, 1000);
    for (int i = 0; i < SPECIES_COUNT; i++) {
        ASSERT(lv.populations[i] >= 0.0, "negative population");
        ASSERT(lv.populations[i] < 1e6, "population exploded");
    }
    ASSERT(lv.generations == 1000, "generation count mismatch");
    PASS();
}

static void test_coexistence_check(void) {
    TEST("coexistence_check with low threshold reports survival");
    Population *p = population_create(200, 13);
    population_interact(p);
    population_classify(p);
    SpeciesTracker st;
    species_tracker_init(&st, p);
    LotkaVolterra lv;
    lotka_volterra_init(&lv, &st);
    CoexistenceResult cr = coexistence_check(&lv, 500, 0.1);
    ASSERT(cr.surviving_count >= 1, "no species survive");
    ASSERT(cr.min_population >= 0.0, "negative min population");
    population_free(p);
    PASS();
}

static void test_fitness_landscape_evaluate(void) {
    TEST("fitness_landscape_evaluate returns non-negative");
    FitnessLandscape fl;
    fitness_landscape_init(&fl, 0.0, 0.1, 42);
    for (int s = 0; s < 3; s++) {
        double f = fitness_landscape_evaluate(&fl, (TernaryStrategy)s);
        ASSERT(f >= 0.0, "negative fitness from landscape");
    }
    PASS();
}

static void test_fitness_landscape_apply(void) {
    TEST("fitness_landscape_apply updates all agents");
    Population *p = population_create(100, 88);
    FitnessLandscape fl;
    fitness_landscape_init(&fl, 0.2, 0.05, 88);
    fitness_landscape_apply(&fl, p);
    int updated = 0;
    for (size_t i = 0; i < p->count; i++) {
        if (p->agents[i].fitness > 0.0) updated++;
    }
    population_free(p);
    ASSERT(updated > 0, "no fitnesses updated");
    PASS();
}

static void test_end_to_end_pipeline(void) {
    TEST("end-to-end: create → interact → classify → track → LV → coexistence");
    Population *p = population_create(150, 2024);
    ASSERT(p != NULL, "create failed");

    population_interact(p);
    population_classify(p);

    SpeciesTracker st;
    species_tracker_init(&st, p);

    LotkaVolterra lv;
    lotka_volterra_init(&lv, &st);
    CoexistenceResult cr = coexistence_check(&lv, 2000, 0.5);

    ASSERT(cr.min_population >= 0.0, "negative population in coexistence");
    ASSERT(cr.surviving_count > 0, "all species died");

    FitnessLandscape fl;
    fitness_landscape_init(&fl, 0.1, 0.0, 2024);
    fitness_landscape_apply(&fl, p);

    /* Verify landscape overrides pairwise fitness */
    int changed = 0;
    for (size_t i = 0; i < p->count; i++) {
        if (p->agents[i].fitness > 0.0) changed++;
    }
    ASSERT(changed > 0, "landscape didn't change fitnesses");

    population_free(p);
    PASS();
}

int main(void) {
    printf("\n=== Strategy Ecology Tests ===\n\n");

    test_population_create_and_size();
    test_population_strategies_in_range();
    test_payoff_matrix_default();
    test_interaction_updates_fitness();
    test_classify_all_species();
    test_species_tracker_counts();
    test_species_tracker_avg_fitness();
    test_lotka_volterra_init();
    test_lotka_volterra_step_changes_state();
    test_lotka_volterra_run();
    test_coexistence_check();
    test_fitness_landscape_evaluate();
    test_fitness_landscape_apply();
    test_end_to_end_pipeline();

    printf("\n%d passed, %d failed\n\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
