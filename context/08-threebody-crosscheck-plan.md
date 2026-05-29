# ThreeBodyDecays Cross-Check Plan

Issue: https://github.com/RUB-EP1/CascadeDecays.jl/issues/3

Reference repositories inspected locally:

- `/Users/mikhailmikhasenko/Documents/DecayModels.CAT/ThreeBodyDecays.jl`
- `/Users/mikhailmikhasenko/Documents/DecayModels.CAT/InstructionalDecayTrees.jl`

## Goal

Build a staged compatibility check between `CascadeDecays.jl` and
`ThreeBodyDecays.jl` for the three-body topology `((1,2),3)`, starting from
aligned kinematics and ending with rotated-frame helicity amplitudes.

The work should produce small convention tests before full amplitude tests. This
is important because failures may come from several independent sources:

- Mandelstam-to-four-vector ordering
- helicity-angle extraction in degenerate or nearly degenerate frames
- particle-2 convention phases
- `ThreeBodyDecays.aligned_amplitude` normalization
- parent orientation phases in non-aligned frames
- Wigner reference-frame rotations

## Current Local State

`CascadeDecays.jl` already has most structural pieces needed for the first pass:

- `DecayTopology(((1,2),3))`
- graph-indexed `CascadeSystem` and `CascadeKinematics`
- `cascade_kinematics(topology, system, objs)` from final-state four-vectors
- generated `helicity_angle_programs(topology)`
- a `ThreeBodyDecays.VertexFunction` adapter in `routed_vertex_amplitude`
- `amplitude(chain, system, x, external_two_λs)` with internal-helicity summation

Current tests pass with:

```sh
julia --project=. test/runtests.jl
```

## Reference Findings

### `ThreeBodyDecays.aligned_four_vectors`

`ThreeBodyDecays.aligned_four_vectors(σs, ms; k)` returns final-state
four-vectors in the parent center-of-mass frame. Particle `k` is aligned with
the negative z axis and all momenta lie in the x-z plane.

For issue A, use `k = 3` to match the topology `((1,2),3)`, where the isobar is
`(1,2)` and particle `3` is the spectator.

### Aligned Amplitude Conventions

`ThreeBodyDecays.aligned_amplitude(dc, σs)` is not only a product of the two
local vertex amplitudes. It also includes:

- `d^J(θ_ij)` for the intermediate spin rotation
- particle-2 convention phase factors at both vertices
- the chain normalization `sqrt(2J + 1)`
- the lineshape and both vertex form factors

This means the first `CascadeDecays.amplitude` comparison should either:

1. include those factors explicitly in the new evaluator conventions, or
2. compare lower-level pieces first and document the known missing factors.

The lower-level route is safer for debugging.

### InstructionalDecayTrees Degeneracy

For aligned `((1,2),3)` kinematics:

- the root decay should have `cosθ = 1` if measuring the `(1,2)` branch after
  boosting to the root helicity frame
- the root azimuth is degenerate when `θ = 0`
- the `(1,2) -> 1,2` decay should have `ϕ = 0` in the aligned x-z plane

The issue suggestion to apply a small positive ZYZ rotation is a good guard
against the undefined root azimuth. A fixed rotation, for example around
`10deg`, should be part of the angle-extraction tests before any amplitude test
depends on `ϕ`.

### Full Orientation Caveat

`ThreeBodyDecays.amplitude(dc, orientation_angles, σs)` applies the parent
orientation through `conj(wignerD_doublearg(two_j0, orientation_angles...))`.

The issue's equation 1 warns that rewriting two D matrices into the full
rotated convention introduces an extra phase involving the spectator helicity
and the inner azimuth. This should be tested after the aligned pieces are
stable, not before.

## Work Schedule

### Stage 0: Baseline and Fixture Selection

Deliverables:

- choose one deterministic mass tuple and one physical `σs`
- choose one chain with simple nonzero recouplings
- store expected raw kinematics in a helper test block

Recommended start:

- use `ThreeBodyMasses(1.1, 2.2, 3.3; m0 = 7.7)` from existing
  `ThreeBodyDecays.jl` tests
- use `x2σs([0.5, 0.3], ms; k = 3)` or another centrally physical point
- fix the random seed only if a generated point is used

### Stage A: Aligned Kinematics Without InstructionalDecayTrees

Purpose: compare the hand-routed `CascadeKinematics` values against
`ThreeBodyDecays` invariant formulas.

Deliverables:

- `test/threebody_compat_tests.jl`
- `CascadeKinematics` built by hand with:
  - root angle `(cosθ = 1, ϕ = 0)`
  - inner angle `(cosθ = cosθij(σs, ms^2; k = 3), ϕ = 0)`
  - line masses squared `[m1^2, m2^2, m3^2, σ3, m0^2]`
- tests that `vertex_masses2`, `vertex_helicities`, and propagator routing
  match the intended `ThreeBodyDecays` chain with `k = 3`

This stage should not depend on four-vectors.

### Stage B: Four-Vector Angle Extraction

Purpose: verify that `helicity_angle_programs` and
`cascade_kinematics(topology, system, objs)` reproduce the Stage A inputs.

Deliverables:

- convert `ThreeBodyDecays.aligned_four_vectors(σs, ms; k = 3)` tuples into
  `FourVectors.FourVector`
- verify invariant masses from `cascade_kinematics`
- verify non-degenerate inner angle against `cosθij`
- document or guard the root azimuth degeneracy
- add a small rotated fixture to verify azimuth behavior away from `θ = 0`

If this stage fails, debug only `helicity_angle_program` and
`InstructionalDecayTrees` path generation before touching amplitude code.

### Stage C: Aligned Amplitude Pieces

Purpose: compare physics factors without reference-frame rotations.

Deliverables:

- compare local `VertexFunction` recouplings for root and inner vertices
- compare form-factor and lineshape calls
- decide where to encode particle-2 convention phases
- decide whether `sqrt(2J + 1)` belongs in `CascadeDecays.amplitude`, a
  `ThreeBodyDecays` compatibility adapter, or documentation as a convention
  difference

Recommended debug order:

1. root vertex only
2. inner vertex only
3. intermediate Wigner small-d only
4. full aligned chain

### Stage D: Rotated Parent Orientation

Purpose: test issue equation 1 without mixing in cross-topology Wigner rotations.

Deliverables:

- generate a non-CMS or explicitly rotated four-vector configuration
- extract parent ZYZ orientation with `InstructionalDecayTrees` tracking
- compare `ThreeBodyDecays.amplitude(dc, orientation_angles, σs)` to the
  corresponding `CascadeDecays` expression plus the documented extra phase
- test both a meson-like spin-0 parent case and a fermion parent case

This stage should remain a single-chain test with the same `k = 3` reference.

### Stage E: Reference-Topology Wigner Rotations

Purpose: generalize beyond the same reference chain.

Deliverables:

- derive graph paths for natural and reference helicity frames
- compare path pairs with `InstructionalDecayTrees.compare_instruction_paths`
- extract `wigner_zyz`, `wigner_zyz_so3`, and/or `wigner_zyz_su2`
- reproduce selected `ThreeBodyDecays` `refζs` cases

This should be the last stage because it depends on every previous convention.

## Suggested Implementation Order

1. Add `test/threebody_compat_tests.jl` and include it from `test/runtests.jl`.
2. Implement Stage A with hand-written `CascadeKinematics`.
3. Implement a small helper that converts `aligned_four_vectors` output to
   `FourVector(px, py, pz; E)`.
4. Implement Stage B and mark root `ϕ` assertions carefully near `θ = 0`.
5. Add local vertex and lineshape comparisons for Stage C.
6. Only then decide whether the main evaluator should change or whether a
   compatibility wrapper is cleaner.
7. Add documentation notes for every convention decision before Stage D/E.

## Early Risks

- The current `routed_vertex_amplitude(::ThreeBodyDecays.VertexFunction, ...)`
  does not visibly include the particle-2 convention phase factors used in
  `ThreeBodyDecays.aligned_amplitude`.
- It also does not include the `sqrt(2J + 1)` chain normalization.
- `helicity_angle_program` currently always uses `ToHelicityFrame`, while
  some helicity paths may need `ToHelicityFrameParticle2` for strict
  Jacob-Wick particle-2 convention compatibility.
- The root azimuth in exactly aligned kinematics is mathematically undefined;
  tests should avoid using it as a hard convention signal.

## Useful Debugging Output

Each stage should print or store, at least temporarily:

- `σs` and masses
- `aligned_four_vectors` converted to `(px, py, pz, E)`
- `CascadeKinematics.line_masses2`
- `CascadeKinematics.vertex_angles`
- root and inner vertex local arguments
- individual phase and normalization factors

These outputs will make issue comments and docs far easier to audit.
