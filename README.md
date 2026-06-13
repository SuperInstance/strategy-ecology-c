# Strategy Ecology — Game Theory Meets Population Dynamics for Ternary Agents

**Strategy ecology** is the study of how strategic agent populations evolve over time when agents interact through game-theoretic payoff structures. This library implements a complete pipeline: agents choose **ternary strategies** (Rock, Paper, Scissors), interact pairwise via a payoff matrix, get classified into **five ecological species** based on fitness, and their population dynamics are modeled by a **Lotka-Volterra system** with Euler integration. The result is a predictive model for whether strategic diversity is stable or will collapse to monoculture.

## Why It Matters

The Rock-Paper-Scissors (RPS) game is the simplest non-trivial game with **no pure-strategy Nash equilibrium** — no single best strategy exists. This makes it the canonical model for **cyclic dynamics**: Rock dominates Scissors, Scissors dominates Paper, Paper dominates Rock, forever. In nature, this exact dynamic appears in side-blotched lizards (orange-blue-yellow throat morphs), in E. coli strains (toxin-producer, toxin-resistant, toxin-sensitive), and in political elections.

In an agent fleet, RPS dynamics arise whenever three strategies are mutually counteracting: a verification agent catches errors from a generation agent, a red-team agent catches biases in the verifier, and the generator exploits gaps in the red-team. None dominates permanently. Understanding whether this tripod is stable — or whether one leg will collapse — is the question this library answers.

The library adds two layers beyond vanilla RPS: (1) **species classification** into five ecological roles (Alpha, Beta, Gamma, Delta, Epsilon), and (2) **Lotka-Volterra population dynamics** on those species, enabling long-horizon coexistence prediction.

## How It Works

### Pipeline Overview

```
1. CREATE     → N agents with random ternary strategies
2. INTERACT   → All-pairs payoff evaluation via RPS matrix
3. CLASSIFY   → Assign each agent to one of 5 species by fitness
4. TRACK      → Count species populations, compute average fitness
5. LV MODEL   → Initialize Lotka-Volterra from species counts
6. CHECK      → Simulate N generations, test coexistence
```

### Step 1: Ternary Strategy Assignment

Each agent carries a strategy from {Rock=0, Paper=1, Scissors=2}. Assignment is uniform random via a seeded LCG (Linear Congruential Generator):

```c
rng_state = seed;
agent[i].strategy = rng_next() % 3;
```

**LCG parameters** (multiplier 1103515245, increment 12345, modulus 2³²) are the ANSI C glibc defaults. Reproducible given the same seed — essential for deterministic testing.

### Step 2: Pairwise Payoff Evaluation

The default RPS payoff matrix:

```
         Rock  Paper  Scissors
Rock     0.5   0.0    1.0
Paper    1.0   0.5    0.0
Scissors 0.0   1.0    0.5
```

A payoff of 1.0 = win, 0.5 = tie, 0.0 = loss. For each pair (i, j), agent i gets `payoff_matrix[s_i][s_j]` and agent j gets `payoff_matrix[s_j][s_i]`.

Fitness is normalized by the number of interactions:

```
fitness_i = (1/(N-1)) · Σ_{j≠i} payoff[si][sj]
```

**Big-O:** `population_interact()` is O(N²) — all-pairs evaluation. For N = 200 agents, this is 19,900 pairwise comparisons. At ~10ns per comparison (array lookup), this is ~200μs — trivially fast on any hardware.

### Step 3: Species Classification

After fitness evaluation, agents are classified into five species based on their relative fitness:

| Species | Condition | Ecological Role |
|---|---|---|
| **Alpha** (0) | Strategy = dominant strategy, top fitness tier | Dominant competitor |
| **Beta** (1) | Strategy = counter to dominant, top fitness tier | Counter-strategist |
| **Gamma** (2) | Strategy = counter-counter, top fitness tier | Cyclic resistor |
| **Delta** (3) | Medium fitness (any strategy) | Generalist / opportunist |
| **Epsilon** (4) | Low fitness (any strategy) | Scavenger / dying |

Classification uses tertiles of the fitness range:

```c
threshold = (max_fit - min_fit) / 3.0;
if (rel_fitness >= 2*threshold) species = STRATEGY_AS_SPECIES;
else if (rel_fitness >= threshold) species = DELTA;
else species = EPSILON;
```

This is **adaptive classification** — the boundaries shift with the population's fitness distribution, not fixed thresholds.

### Step 4: Lotka-Volterra Dynamics

Species populations feed into a 5-species Lotka-Volterra system:

```
dNᵢ/dt = rᵢ · Nᵢ · (1 - Nᵢ/Kᵢ) + Σ_{j≠i} Aᵢⱼ · Nᵢ · Nⱼ - δ · Nᵢ
```

Three terms:
1. **Logistic growth:** `rᵢ · Nᵢ · (1 - Nᵢ/Kᵢ)` — exponential growth bounded by carrying capacity
2. **Interaction:** `Σⱼ Aᵢⱼ · Nᵢ · Nⱼ` — species j affects species i (positive = mutualism, negative = competition)
3. **Decay:** `δ · Nᵢ` — universal death rate

The interaction matrix encodes **cyclic dominance** among Alpha/Beta/Gamma:

```
          ALP     BET     GAM
ALP    {-0.01,  0.005, -0.008}   ← Alpha helped by Beta, hurt by Gamma
BET    {-0.008, -0.01,   0.005}  ← Beta helped by Gamma, hurt by Alpha
GAM    { 0.005, -0.008, -0.01}   ← Gamma helped by Alpha, hurt by Beta
```

This cyclic structure (A→B→C→A) is what produces stable oscillations rather than competitive exclusion.

**Integration:** Forward Euler with dt = 0.01. Per step: O(S²) where S = number of species (5), so 25 multiply-adds per step — sub-microsecond.

### Step 5: Coexistence Check

After simulating N generations, the library checks whether all species survive above a threshold:

```c
typedef struct {
    bool all_survive;
    double min_population;
    int surviving_count;
    double final_populations[SPECIES_COUNT];
} CoexistenceResult;
```

If `all_survive = false`, the system is heading toward **competitive exclusion** — one or more species will go extinct. If `surviving_count < SPECIES_COUNT`, partial collapse has occurred.

### Fitness Landscape

The `FitnessLandscape` provides an alternative fitness evaluation — instead of pairwise RPS interaction, fitness is computed from a Gaussian landscape:

```
fitness(s) = peak_s · exp(-(shifted_s)² / (2 · width_s²)) + noise
```

Where `shifted_s = peak_s - shift · (s - 1)`. The `shift` parameter models **environmental change**: when shift = 0, all strategies have equal fitness peaks. When shift = 1, strategies are offset, modeling a changing environment where different strategies are favored at different times.

## Quick Start

```bash
# Build
gcc -o test_ecology tests/test_ecology.c src/strategy_ecology.c -lm -Wall -O2

# Run all 14 tests
./test_ecology
```

```c
#include "strategy_ecology.h"

int main(void) {
    // Create 200 agents with seed 42
    Population *pop = population_create(200, 42);

    // Evaluate all pairwise interactions
    population_interact(pop);

    // Classify into species
    population_classify(pop);

    // Track species counts
    SpeciesTracker st;
    species_tracker_init(&st, pop);
    for (int i = 0; i < SPECIES_COUNT; i++) {
        printf("Species %d: %d agents, avg fitness %.3f\n",
               i, species_tracker_count(&st, (SpeciesID)i),
               species_tracker_avg_fitness(&st, (SpeciesID)i));
    }

    // Initialize Lotka-Volterra from species counts
    LotkaVolterra lv;
    lotka_volterra_init(&lv, &st);

    // Simulate 2000 generations
    CoexistenceResult cr = coexistence_check(&lv, 2000, 0.5);
    printf("All survive: %s, surviving: %d/%d\n",
           cr.all_survive ? "yes" : "no",
           cr.surviving_count, SPECIES_COUNT);

    // Apply fitness landscape with environmental shift
    FitnessLandscape fl;
    fitness_landscape_init(&fl, 0.2, 0.05, 42);
    fitness_landscape_apply(&fl, pop);

    population_free(pop);
    return 0;
}
```

## API

### Types

```c
typedef enum { STRAT_ROCK=0, STRAT_PAPER=1, STRAT_SCISSORS=2 } TernaryStrategy;

typedef enum {
    SPECIES_ALPHA=0,   // dominant strategy, high fitness
    SPECIES_BETA=1,    // counter-strategy, high fitness
    SPECIES_GAMMA=2,   // counter-counter, high fitness
    SPECIES_DELTA=3,   // generalist, medium fitness
    SPECIES_EPSILON=4, // scavenger, low fitness
    SPECIES_COUNT=5
} SpeciesID;
```

### Population

| Function | Description | Complexity |
|---|---|---|
| `population_create(n, seed)` | Allocate N agents with random strategies | O(N) |
| `population_interact(pop)` | All-pairs payoff evaluation | O(N²) |
| `population_classify(pop)` | Assign species by fitness tertile | O(N) |
| `population_free(pop)` | Deallocate | O(1) |

### SpeciesTracker

| Function | Description |
|---|---|
| `species_tracker_init(st, pop)` | Count species from classified population |
| `species_tracker_count(st, id)` | Population count for species `id` |
| `species_tracker_avg_fitness(st, id)` | Mean fitness for species `id` |

### LotkaVolterra

| Function | Description | Complexity |
|---|---|---|
| `lotka_volterra_init(lv, st)` | Initialize from species tracker | O(S²) |
| `lotka_volterra_step(lv)` | One Euler integration step | O(S²) |
| `lotka_volterra_run(lv, n)` | Run N steps | O(N · S²) |

### Coexistence & Fitness

| Function | Description |
|---|---|
| `coexistence_check(lv, generations, threshold)` | Simulate and check survival |
| `fitness_landscape_init(fl, shift, noise, seed)` | Configure Gaussian landscape |
| `fitness_landscape_evaluate(fl, strategy)` | Fitness for one strategy |
| `fitness_landscape_apply(fl, pop)` | Override population fitness from landscape |

## Architecture Notes

This library is the **on-device strategy engine** for the SuperInstance constellation. In the conservation law **γ + η = C**, strategy ecology predicts how γ (generation energy) should be allocated across competing agent strategies. When three strategies coexist, γ is diversified — the system is resilient to any single strategy's failure. When one strategy dominates, γ is concentrated — the system is efficient but fragile.

The library pairs with `lotka-volterra-agents-c` (the general dynamics engine) to form the complete ecological prediction layer for edge deployment. See the [SuperInstance Architecture](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## References

1. Sinervo, B. & Lively, C. "The Rock-Paper-Scissors Game and the Evolution of Alternative Male Strategies" (1996) — empirical RPS dynamics in lizards, *Nature* 380:240–243
2. May, R. M. *Stability and Complexity in Model Ecosystems* (1973) — Lotka-Volterra stability analysis
3. Hofbauer, J. & Sigmund, K. *Evolutionary Games and Population Dynamics* (1998) — game theory ↔ dynamical systems correspondence

## License

MIT
