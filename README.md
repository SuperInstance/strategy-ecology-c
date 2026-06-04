# strategy-ecology-c

C implementation of strategy ecology — Lotka-Volterra style species dynamics for ternary agents.

## Components

- **Population**: N agents with ternary strategies (Rock/Paper/Scissors)
- **Interaction**: pairwise payoff matrix evaluation
- **SpeciesTracker**: count agents in each of 5 species categories
- **LotkaVolterra**: Euler integration of species dynamics (growth, interaction, decay)
- **CoexistenceCheck**: verify all species survive after N generations
- **FitnessLandscape**: evaluate fitness in configurable landscape

## Build & Test

```bash
gcc -o test_ecology tests/test_ecology.c src/strategy_ecology.c -lm -Wall -O2
./test_ecology
```

## License

MIT
