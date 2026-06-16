# CascadeDecays.jl

`CascadeDecays.jl` provides topology-aware routing for sequential decays: it stores a flat binary cascade graph, resolves bracket-addressed propagators and vertices, and maps event-dependent masses and helicity angles into local amplitude building blocks.

## Getting started

- Define a cascade with `DecayTopology(...)`.
- Attach propagators and vertices with `DecayChain(...)`.
- Build runtime kinematics with `CascadeKinematics(topology, objs)`.
- Evaluate amplitudes with `amplitude(chain, system, x)` (full external-helicity array) or `amplitude(..., external_two_λs)` (one component).

The main user-facing walkthrough is the tutorial page generated from the Quarto source in `docs/integration_4body_b2ddKpi.qmd`.

For a practical multi-channel model construction sketch, see [Lb to Lc 3pi model](@ref lb2lc3pi_model).

For a three-body cross-check against `ThreeBodyDecays.jl`, including reference-topology Wigner rotations and overall sign conventions, see [CascadeDecays vs DPD](@ref cascade_vs_dpd).

For line ids, vertex indices, and the incidence matrix, see [Internal notation](@ref notation).
