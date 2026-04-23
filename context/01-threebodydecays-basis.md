# Basis From `ThreeBodyDecays.jl`

## Scope of the review

The main source inspected was the upstream repository:

- `src/decay_channel.jl`
- `src/recouplings.jl`
- `src/wigner_d_matrix.jl`
- `src/wigner_rotations.jl`
- `src/mass_types.jl`
- `src/spin_types.jl`
- `src/system_types.jl`
- `docs/src/10-energy-dependence.md`
- `docs/src/11-quantization-reference.md`

The important conclusion is that `ThreeBodyDecays.jl` already separates the physics into three layers:

1. kinematic system data
2. local building blocks for a two-body decay step
3. an evaluation recipe for a complete chain

That separation should be preserved in `CascadeDecays.jl`, but generalized from the fixed three-body setting to arbitrary cascade graphs.

## What the current package really stores

### 1. External system data

`ThreeBodySystem` stores only the external masses and spins:

- `ms`: a `MassTuple` with `(m1, m2, m3, m0)`
- `two_js`: a `SpinTuple` with `(two_h1, two_h2, two_h3, two_h0)`

This is a strong signal for the new package:

- topology and quantum-number payload should be stored separately
- event-dependent internal kinematics should not be baked into the type

### 2. Vertex content

The local vertex object is `VertexFunction{R,F}` with:

- `h::R`, where `R <: Recoupling`
- `ff::F`, a callable form factor

The important part is that the vertex is not only a coupling rule. It is:

- spin/helicity recoupling
- optional energy dependence through a form factor

Concrete recoupling variants already present:

- `NoRecoupling`
- `ParityRecoupling`
- `RecouplingLS`

This means the new package should likely distinguish:

- a vertex payload type for quantum-number coupling
- an optional vertex-attached dynamical factor

### 3. Propagator / lineshape content

`ThreeBodyDecays.jl` does not define a dedicated abstract `Propagator` type in the inspected code.
Instead, the propagating intermediate state is represented by a callable object stored as `Xlineshape`.

From the docs:

- production form factor is attached to `HRk`
- propagator is a callable of one variable, `σk`
- decay form factor is attached to `Hij`

So the physical factorization is:

`vertex * propagator * vertex`

not

`propagator hidden inside a vertex`

That matches your desired design very well. In the new package it is worth introducing an explicit abstract interface such as:

- `AbstractLineshape`
- `AbstractVertex`

even though the upstream package does not yet formalize them that way.

### 4. Chain container

The current chain object is:

- `DecayChain`

with fields:

- `k`: topology label for the spectator choice
- `two_j`: spin of the intermediate isobar
- `Xlineshape`: propagator / lineshape
- `HRk`: parent-to-isobar vertex
- `Hij`: isobar-to-daughters vertex
- `tbs`: full external system

This is a single fixed three-body topology container. It already shows the pattern we want to generalize:

- one propagating line
- two local vertices
- a fixed connectivity rule derived from the chain type itself

## What should be inherited conceptually

### Separation of responsibilities

The upstream package cleanly separates:

- masses/spins of external particles
- vertex recouplings
- form factors
- propagator
- Wigner rotations
- topology-specific evaluation logic

For `CascadeDecays.jl`, the same separation should stay, but the topology logic should move from hard-coded chain formulas to a flat graph container plus derived evaluation order.

## Wigner machinery

There are two distinct Wigner ingredients in `ThreeBodyDecays.jl`:

### Wigner small-d / D functions

`wignerD_doublearg` is imported from `PartialWaveFunctions`.
`wigner_d_matrix.jl` mostly provides convenience broadcasting wrappers.

So the real design lesson is:

- the package depends on an external implementation of the special functions
- the package-level job is to organize the indices and conventions, not necessarily to reimplement special functions from scratch

This suggests `CascadeDecays.jl` should initially wrap an existing reliable Wigner implementation and focus on:

- argument conventions
- phase conventions
- index ordering
- caching / evaluation strategy

### Wigner rotation bookkeeping

`wigner_rotations.jl` is not just special functions. It encodes how different chain-dependent quantization frames are related.

This is conceptually important for the new package:

- graph topology alone is not enough
- we also need a notion of ordered boosts / decay-frame transitions
- those transitions are what determine the Wigner rotations appearing between local amplitudes

For a general cascade package, this likely becomes a derived structure from the graph:

- path from mother to each external leg
- path comparisons between alternative coupling schemes
- local rotation objects computed from path mismatch

## Concrete lessons for the new package

### 1. Introduce explicit abstract interfaces

The upstream code uses plain callables in several places. For a new general package, explicit interfaces will help:

- `AbstractLineshape`
- `AbstractVertex`
- possibly `AbstractFormFactor`
- possibly `AbstractRecoupling`

### 2. Do not encode topology recursively in nested Julia objects

In `ThreeBodyDecays.jl`, topology is still small enough to be implicit in the chain definition. For arbitrary cascades that will not scale.

Your flat graph representation is the right generalization.

### 3. Keep dynamics attached to the correct objects

- vertices carry recouplings and maybe form factors
- propagators belong to internal lines only
- external legs should not require propagators

### 4. Keep external system parameters separate from event-dependent internals

In the upstream package the invariant `σk` is runtime kinematics, not stored in the chain object.
That same idea should be preserved for every internal line in a general cascade.

## Recommended translation into `CascadeDecays.jl`

The direct conceptual translation is:

- `ThreeBodySystem` becomes a more general external-state descriptor
- `DecayChain` becomes a flat graph topology container
- `VertexFunction` informs the future `AbstractVertex` payload
- `Xlineshape` informs the future `AbstractLineshape` payload
- Wigner wrappers remain utility-level physics functions
- frame-transport and computational ordering should be derived from graph paths

## One especially important finding

The upstream package does not truly treat the intermediate mass as a stored constant.
It stores the external masses, while the intermediate line is evaluated at the kinematic invariant `σk`.

For the new package this is crucial:

- intermediate lines are topological objects with quantum numbers and dynamics
- their instantaneous invariant masses belong to phase-space points or derived event kinematics
- they should not be frozen as numeric fields in the canonical topology container
