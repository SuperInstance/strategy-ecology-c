# Future Integration: strategy-ecology-c

## Current State
C implementation of strategy ecology — Lotka-Volterra style species dynamics for ternary agents. Models how different ternary strategies compete, cooperate, and coexist in an ecosystem.

## Integration Opportunities

### With lotka-volterra-agents-c
Shared mathematical foundation. `lotka-volterra-agents-c` is the general engine; `strategy-ecology-c` is the strategy-specific application. Together they form a complete ecological simulation for ternary agents on embedded devices.

### With ternary-game-theory
Strategic interactions in ecology ARE games. The payoff matrix from `ternary-game-theory` defines the interaction coefficients for the Lotka-Volterra dynamics in `strategy-ecology-c`. Nash equilibrium of the game corresponds to the ecological equilibrium.

### With ternary-cell
Cell types with different strategies form an ecology. The `vibe` phase updates cell states based on neighbor interactions — this IS strategy ecology. The C port runs the ecological dynamics on-device for real-time cell population management.

## Potential in Mature Systems
In room-as-codespace, rooms host strategy ecologies. Different ensign types compete for the right to handle tasks. The C port predicts which strategy will dominate given current conditions, enabling the room to pre-load the winning ensign before the competition resolves.

## Cross-Pollination Ideas
- Strategy coexistence as room diversity metric — healthy rooms have multiple coexisting strategies
- Ecological equilibrium as room stability — rooms at ecological equilibrium are predictable
- Invasive strategy detection — a new strategy displacing natives indicates room configuration drift

## Dependencies for Next Steps
- Integration with lotka-volterra-agents-c for combined simulation
- FFI bindings for Rust interop
- Calibration data from ternary-cell population observations
