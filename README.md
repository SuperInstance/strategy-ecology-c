# strategy-ecology-c

C implementation of strategy species ecology for ternary agents.

Five species: **Explorer**, **Diplomat**, **Marksman**, **Climber**, **Prospector**.

## Components

- **Species enum & traits** — typical win rates and entropy per species
- **Population** — species counts, Shannon diversity, Simpson index
- **LotkaVolterra** — Euler simulation with interaction matrix
- **EcologicalResilience** — species survival and resilience index
- **SpeciesClassify** — classify agents into species by behavior profile

## Build & Test

```bash
gcc -o test_ecology tests/test_ecology.c src/strategy_ecology.c -lm -Wall -O2
./test_ecology
```

## License

MIT
