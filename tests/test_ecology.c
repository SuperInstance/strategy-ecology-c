#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../src/strategy_ecology.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT_EQ_DBL(a, b, eps) do { \
    double _a = (a), _b = (b), _eps = (eps); \
    if (fabs(_a - _b) > _eps) { \
        printf("  FAIL %s:%d: %.6f != %.6f (eps=%.6f)\n", __FILE__, __LINE__, _a, _b, _eps); \
        tests_failed++; \
    } else { tests_passed++; } \
} while(0)

#define ASSERT_EQ_INT(a, b) do { \
    int _a = (a), _b = (b); \
    if (_a != _b) { \
        printf("  FAIL %s:%d: %d != %d\n", __FILE__, __LINE__, _a, _b); \
        tests_failed++; \
    } else { tests_passed++; } \
} while(0)

#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { \
        printf("  FAIL %s:%d: assertion false\n", __FILE__, __LINE__); \
        tests_failed++; \
    } else { tests_passed++; } \
} while(0)

#define TEST(name) printf("── %s ──\n", name)

/* ── Tests ────────────────────────────────────────────────────── */

static void test_species_traits_not_null(void) {
    TEST("species_traits returns non-null for valid species");
    for (int i = 0; i < SPECIES_COUNT; i++) {
        const SpeciesTraits *t = species_traits((Species)i);
        ASSERT_TRUE(t != NULL);
        ASSERT_TRUE(t->win_rate > 0.0 && t->win_rate < 1.0);
        ASSERT_TRUE(t->entropy > 0.0);
        ASSERT_TRUE(t->name != NULL);
    }
}

static void test_species_traits_invalid(void) {
    TEST("species_traits returns NULL for invalid species");
    ASSERT_TRUE(species_traits((Species)-1) == NULL);
    ASSERT_TRUE(species_traits((Species)SPECIES_COUNT) == NULL);
}

static void test_population_total(void) {
    TEST("population_total");
    unsigned int counts[SPECIES_COUNT] = {10, 20, 30, 15, 25};
    Population p = population_create(counts);
    ASSERT_EQ_INT(population_total(&p), 100);
}

static void test_population_total_zero(void) {
    TEST("population_total zero");
    unsigned int counts[SPECIES_COUNT] = {0, 0, 0, 0, 0};
    Population p = population_create(counts);
    ASSERT_EQ_INT(population_total(&p), 0);
}

static void test_population_fraction(void) {
    TEST("population_fraction");
    unsigned int counts[SPECIES_COUNT] = {20, 20, 20, 20, 20};
    Population p = population_create(counts);
    ASSERT_EQ_DBL(population_fraction(&p, SPECIES_EXPLORER), 0.2, 1e-9);
    ASSERT_EQ_DBL(population_fraction(&p, SPECIES_MARKSMAN), 0.2, 1e-9);
}

static void test_population_shannon_uniform(void) {
    TEST("population_shannon uniform → log(5)");
    unsigned int counts[SPECIES_COUNT] = {100, 100, 100, 100, 100};
    Population p = population_create(counts);
    double expected = log(5.0);
    ASSERT_EQ_DBL(population_shannon(&p), expected, 1e-9);
}

static void test_population_shannon_zero(void) {
    TEST("population_shannon empty → 0");
    unsigned int counts[SPECIES_COUNT] = {0, 0, 0, 0, 0};
    Population p = population_create(counts);
    ASSERT_EQ_DBL(population_shannon(&p), 0.0, 1e-9);
}

static void test_population_simpson_uniform(void) {
    TEST("population_simpson uniform");
    unsigned int counts[SPECIES_COUNT] = {100, 100, 100, 100, 100};
    Population p = population_create(counts);
    double expected = 1.0 - 5.0 * 0.04; /* 1 - 5*(1/5)² = 0.8 */
    ASSERT_EQ_DBL(population_simpson(&p), expected, 1e-9);
}

static void test_population_simpson_dominant(void) {
    TEST("population_simpson single dominant species");
    unsigned int counts[SPECIES_COUNT] = {96, 1, 1, 1, 1};
    Population p = population_create(counts);
    double s = population_simpson(&p);
    ASSERT_TRUE(s < 0.1); /* very low diversity */
}

static void test_lotka_volterra_step_positive(void) {
    TEST("lotka_volterra_step keeps small populations positive");
    unsigned int counts[SPECIES_COUNT] = {10, 10, 10, 10, 10};
    Population p = population_create(counts);
    LotkaVolterra lv = lotka_volterra_default(0.1);
    Population next = lotka_volterra_step(&lv, &p);
    for (int i = 0; i < SPECIES_COUNT; i++) {
        ASSERT_TRUE(next.counts[i] > 0);
    }
}

static void test_lotka_volterra_simulate_many_steps(void) {
    TEST("lotka_volterra_simulate 1000 steps doesn't collapse to zero");
    unsigned int counts[SPECIES_COUNT] = {50, 50, 50, 50, 50};
    Population p = population_create(counts);
    LotkaVolterra lv = lotka_volterra_default(0.05);
    Population result = lotka_volterra_simulate(&lv, &p, 1000);
    ASSERT_TRUE(population_total(&result) > 0);
}

static void test_lotka_volterra_equilibrium(void) {
    TEST("lotka_volterra large populations decline (competition)");
    unsigned int counts[SPECIES_COUNT] = {1000, 1000, 1000, 1000, 1000};
    Population p = population_create(counts);
    LotkaVolterra lv = lotka_volterra_default(0.1);
    Population next = lotka_volterra_step(&lv, &p);
    /* With high counts, interaction term > 1 → dn/dt negative → population should shrink */
    ASSERT_TRUE(next.counts[0] < p.counts[0]);
}

static void test_ecological_resilience_full(void) {
    TEST("ecological_resilience all survive");
    unsigned int counts[SPECIES_COUNT] = {20, 20, 20, 20, 20};
    Population p = population_create(counts);
    EcologicalResilience er = ecological_resilience(&p, &p);
    ASSERT_EQ_INT(er.surviving_count, 5);
    ASSERT_TRUE(er.resilience_index > 0.9);
}

static void test_ecological_resilience_partial(void) {
    TEST("ecological_resilience some extinct");
    unsigned int before[SPECIES_COUNT] = {20, 20, 20, 20, 20};
    unsigned int after[SPECIES_COUNT]  = {10, 0, 15, 0, 5};
    Population b = population_create(before);
    Population a = population_create(after);
    EcologicalResilience er = ecological_resilience(&b, &a);
    ASSERT_EQ_INT(er.surviving_count, 3);
    ASSERT_EQ_INT(er.surviving[SPECIES_DIPLOMAT], 0);
    ASSERT_EQ_INT(er.surviving[SPECIES_CLIMBER], 0);
    ASSERT_TRUE(er.resilience_index < 1.0);
}

static void test_species_classify_marksman(void) {
    TEST("species_classify high win rate → Marksman");
    AgentBehavior b = { .win_rate = 0.68, .entropy = 0.7, .cooperation_rate = 0.1, .risk_tolerance = 0.5 };
    ASSERT_EQ_INT(species_classify(&b), SPECIES_MARKSMAN);
}

static void test_species_classify_diplomat(void) {
    TEST("species_classify high cooperation → Diplomat");
    AgentBehavior b = { .win_rate = 0.50, .entropy = 1.4, .cooperation_rate = 0.9, .risk_tolerance = 0.1 };
    ASSERT_EQ_INT(species_classify(&b), SPECIES_DIPLOMAT);
}

static void test_species_classify_climber(void) {
    TEST("species_classify high risk → Climber");
    AgentBehavior b = { .win_rate = 0.38, .entropy = 2.3, .cooperation_rate = 0.2, .risk_tolerance = 0.95 };
    ASSERT_EQ_INT(species_classify(&b), SPECIES_CLIMBER);
}

static void test_lotka_volterra_default_dt(void) {
    TEST("lotka_volterra_default stores dt");
    LotkaVolterra lv = lotka_volterra_default(0.42);
    ASSERT_EQ_DBL(lv.dt, 0.42, 1e-12);
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void) {
    test_species_traits_not_null();
    test_species_traits_invalid();
    test_population_total();
    test_population_total_zero();
    test_population_fraction();
    test_population_shannon_uniform();
    test_population_shannon_zero();
    test_population_simpson_uniform();
    test_population_simpson_dominant();
    test_lotka_volterra_step_positive();
    test_lotka_volterra_simulate_many_steps();
    test_lotka_volterra_equilibrium();
    test_ecological_resilience_full();
    test_ecological_resilience_partial();
    test_species_classify_marksman();
    test_species_classify_diplomat();
    test_species_classify_climber();
    test_lotka_volterra_default_dt();

    printf("\n══ Results: %d passed, %d failed ══\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
