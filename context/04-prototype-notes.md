# Prototype Notes

## What the prototype already gets right

The prototype expresses several good design instincts:

- explicit two-body local objects
- separation of masses and spins
- separation of the local decay object from the full cascade
- explicit use of Wigner `D` functions
- explicit summation over internal helicities

This is the right local physics decomposition.

## Strong parts worth preserving

### `TwoBodyMasses`

This is a useful local kinematic descriptor for a single decay vertex:

- parent mass `m0`
- daughter masses `m1`, `m2`

It is a good local evaluation object, but it should probably be runtime kinematics rather than part of the immutable canonical topology container.

### `TwoBodySpins`

This is a good idea.
The use of doubled spins/helicities matches the conventions in `ThreeBodyDecays.jl` and avoids half-integer instability in indexing and dispatch.

### `TwoBodyDecay`

This is also a good local object:

- a two-body system
- a vertex function

Conceptually it is close to one local decay step in the cascade.

### `SimpleCascade`

This is a useful prototype for evaluation logic:

- multiply local amplitudes
- sum over internal helicities

That is exactly the algebra the general graph evaluator will need to automate.

## Where the prototype is still too local

The prototype hard-codes:

- a chain length of three local decays
- explicit intermediate helicity names
- explicit ordered tuple input `(Ω1, Ω2, Ω3)`
- explicit multiplication order in handwritten code

That is perfect for testing ideas, but it must become derived from the graph container.

## Mapping prototype concepts to the target package

Suggested translation:

- `TwoBodyMasses` -> local runtime kinematics of one vertex evaluation
- `TwoBodySpins` -> reusable quantum-number helper type
- `TwoBodySystem` -> local edge/vertex evaluation context
- `TwoBodyDecay` -> candidate basis for a concrete `AbstractVertex` implementation
- `SimpleCascade` -> replaced by generic graph traversal + contraction engine

## Important observation about `amplitude`

In the prototype, the local amplitude contains:

- local recoupling
- local form factor
- local Wigner `D`

That is a sensible starting point, but in the full package we should keep the factorization explicit:

- vertex coupling
- vertex form factor
- propagator on internal line
- Wigner rotation / decay-angle factor

This keeps the design closer to the upstream package and makes substitution easier.

## Recommended next prototype step

Before implementing a full generic evaluator, a good intermediate target would be:

1. define a flat topology for a single four-body example
2. write a function that reconstructs all local two-body decay steps from the graph
3. evaluate the cascade by iterating over those steps
4. only then generalize the contraction over internal helicities

That would give an incremental path from the current handwritten `SimpleCascade` to the final package architecture.
