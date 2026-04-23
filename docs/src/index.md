# CascadeDecays.jl

`CascadeDecays.jl` provides topology-aware routing for sequential decays: it stores a flat binary cascade graph, resolves bracket-addressed propagators and vertices, and maps event-dependent masses and helicity angles into local amplitude building blocks.

## Getting started

- Define a cascade with `DecayTopology(...)`.
- Attach propagators and vertices with `DecayChain(...)`.
- Build runtime kinematics with `cascade_kinematics(...)`.
- Evaluate amplitudes with `amplitude(...)`.

The main user-facing walkthrough is the tutorial page generated from the literate source in `docs/integration_4body_b2ddKpi.jl`.
