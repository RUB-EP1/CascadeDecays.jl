# [Amplitude computation](@id amplitude_computation)

```@meta
CurrentModule = CascadeDecays
EditURL = "../src/amplitude-computation.md"
```

This page describes how [`amplitude`](@ref) assembles a decay-chain amplitude from
the static [`DecayChain`](@ref) payload and runtime [`DecayChainKinematics`](@ref).
It complements [Topology and numbering](@ref notation) (indices and axis layout)
and [Routing four-vectors](@ref kinematic_tasks) (how local angles and masses are
obtained). Tutorials show end-to-end model building; this page is the evaluation
reference.

## Overview

For a single [`DecayChain`](@ref) and one phase-space point `x`, evaluation follows
the same aligned-style chain used in `ThreeBodyDecays.aligned_amplitude`:

1. Build a **vertex factor** `V_v` at every topology vertex `v` (local Wigner
   rotation × coupling × optional form factor).
2. Multiply vertex factors into a **line-indexed buffer**
   [`line_amplitude_tensor`](@ref) `F` (one helicity axis per line).
3. **Sum internal helicities** on propagating lines with
   [`external_helicity_amplitude`](@ref).
4. Multiply by the **propagator product** and a **spin normalisation** factor.

The public entry point is

```julia
A = amplitude(chain, system, x)
```

which returns a complex array with one axis per external helicity: final-state lines
in [`final_line_inds`](@ref) order, then the root helicity last. A single component
is selected with `amplitude(chain, system, x, external_two_λs)`.

!!! note "Aligned vs helicity-frame output"
    `amplitude(chain, system, x)` does **not** apply the four final-state Wigner
    `d`-rotations that `ThreeBodyDecays.amplitude(dc, σs)` adds after the aligned
    chain. The result is an **aligned** amplitude in the natural frames of this
    topology. To compare or add channels with different topologies coherently, use
    [`KinematicTask`](@ref) / [`KinematicPoint`](@ref) and the Wigner-alignment paths
    documented in [Routing four-vectors](@ref kinematic_tasks). See also
    [Cross-checking with ThreeBodyDecays](@ref cascade_vs_dpd).

Schematically,

```math
A_{\lambda_{\mathrm{ext}}} =
\sqrt{\prod_R (2J_R+1)} \;
\prod_{R} X_R(\sigma_R) \;
\sum_{\lambda_{\mathrm{int}}}
\prod_v V_v^{\lambda_0 \lambda_1 \lambda_2},
```

where `R` runs over internal (propagating) lines, `σ_R` is the invariant mass
squared on line `R`, and the sum runs over the helicities carried by those lines.

## Propagators

Each **internal line** carries one propagator payload in `chain.propagators`,
paired with the line id from [`propagating_line_inds`](@ref). Final-state lines
and the root line have no propagator.

A [`Propagator`](@ref) stores

- `two_j` — twice the spin of the resonance on that line (or `SpinParity` metadata
  when built from LS helpers);
- `lineshape` — a callable model (typically [`AbstractLineshape`](@ref)), evaluated
  at the **event-dependent** invariant mass squared `σ` on that line.

At evaluation time the propagator contribution is

```math
X_R(\sigma_R) = \mathcal{X}_R(\sigma_R),
```

i.e. the lineshape callable applied to `line_invariant(x, line_ind)`. The product
over all propagating lines is `routed_propagator_product(chain, x)`.

Additionally, each propagating line contributes a factor `\sqrt{2J_R+1}` in the
overall normalisation (matching `ThreeBodyDecays.aligned_amplitude`).

## Vertices

Each binary vertex `v` connects a parent line `0` and two children `1` and `2`
([`vertex_line_inds`](@ref)). The local payload is a [`Vertex`](@ref): a recoupling
scheme `h` (e.g. `RecouplingLS` from `ThreeBodyDecays`) and an optional
form factor `ff(m_0^2, m_1^2, m_2^2)`.

### Helicity coupling `H` vs recoupling `h`

The helicity coupling used in the amplitude is **not** the raw recoupling object
`h`. Following the Jacob–Wick **particle-2 convention** (as in
`ThreeBodyDecays.jl`), the package defines

```math
H_{\lambda_1 \lambda_2} =
h_{\lambda_1 \lambda_2} \;
\mathrm{phase}(j_2 - \lambda_2) \;
\mathrm{ff}(m_0^2, m_1^2, m_2^2),
```

where `j_2` is the spin of **child 2** (the right child in the bracket tree). The
phase factor is

```math
\mathrm{phase}(j_2 - \lambda_2) = (-1)^{j_2 - \lambda_2}.
```

In the implementation (`_particle_two_phase` in `evaluation.jl`), twice-helicity
integers are used:

```math
\mathrm{phase}(2j_2 - 2\lambda_2) =
\begin{cases}
(-1)^{(2j_2 - 2\lambda_2)/2} & 2j_2 - 2\lambda_2 \text{ even}, \\
\text{error} & \text{otherwise},
\end{cases}
```

which is equivalent to `(-1)^{j_2 - \lambda_2}` whenever `j_2 - \lambda_2` is an
integer (as for physical helicity states).

**Why split `h` and the phase?** The recoupling amplitude `h` has the parity
properties required for LS coupling: [`minimal_ls_decay_chain`](@ref) and
[`possible_vertex_couplings`](@ref) build vertices from allowed `(L, S)` pairs
attached to `h`. The particle-2 phase is a fixed sign factor that implements the
Jacob–Wick choice for quantising child 2; it is applied automatically at evaluation
and should **not** be folded into user-defined `h` payloads.

The kinematic side of the same convention uses `ToHelicityFrameParticle2` (from
`InstructionalDecayTrees`) on Wigner-alignment paths when the transported line is
`child2` at a vertex (see [Routing four-vectors](@ref kinematic_tasks)).

### Full vertex factor

For each helicity triplet `(λ_0, λ_1, λ_2)` on the three lines at vertex `v`, the
vertex factor is

```math
V_v^{\lambda_0 \lambda_1 \lambda_2} =
D^{j_0*}_{\lambda_0,\, \lambda_1 - \lambda_2}(\phi, \theta, 0)\;
H_{\lambda_1 \lambda_2},
```

where `(ϕ, cosθ)` are the local decay angles at that vertex (`vertex_angles(x, vertex_ind)`),
and the Wigner `D`-functions are evaluated using `PartialWaveFunctions.jl`. Entries
with `|\lambda_1 - \lambda_2| > j_0` are zero.

Vertex factors are built as dense arrays (analogous to `VRk` / `Vij` in
`ThreeBodyDecays.aligned_amplitude`) and multiplied into the line-indexed tensor
[`line_amplitude_tensor`](@ref).

### Routed single-vertex helper

For tests and diagnostics, `routed_vertex_amplitude` evaluates one vertex with
masses, helicities, spins, and angles routed from a full cascade context:

```julia
routed_vertex_amplitude(chain, system, x, two_λs, vertex_ind)
```

This is the same local rotation × `H` logic without the chain product or propagators.

## Related API

| Symbol | Role |
|--------|------|
| [`amplitude`](@ref) | Full external-helicity array for one chain and kinematics |
| [`line_amplitude_tensor`](@ref) | Line-indexed product of vertex factors (before internal sum) |
| [`external_helicity_amplitude`](@ref) | Sum internal propagator helicities |
| `routed_propagator_product` | Product of lineshape values on internal lines |
| `routed_vertex_amplitude` | Single-vertex amplitude with routing |
| [`minimal_ls_decay_chain`](@ref) | Build a chain with LS recoupling on `h` |

## Further reading

- [Topology and numbering](@ref notation) — bracket child order fixes `child1` /
  `child2` and therefore which line receives the particle-2 phase.
- [Cross-checking with ThreeBodyDecays](@ref cascade_vs_dpd) — numerical comparison
  of aligned and helicity-frame conventions for a three-body reference.
- `ThreeBodyDecays.jl` — [Quantization reference frames](https://rub-ep1.github.io/ThreeBodyDecays.jl/dev/quantization_reference/) and `aligned_amplitude` in `decay_channel.jl`.
