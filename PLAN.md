# PLAN

## Purpose

`CascadeDecays.jl` should become a Julia package for general hadronic cascade-decay computations built around:

- typed local physics payloads
- a flat line-vertex incidence graph
- derived evaluation order and frame transport

The package should inherit the successful local factorization from `ThreeBodyDecays.jl` while removing the fixed three-body topology assumption.

## Design principles

1. Keep topology flat, not recursively nested in Julia types.
2. Separate topology, dynamics, and runtime kinematics.
3. Attach propagators to internal lines only.
4. Use doubled-spin conventions for quantum numbers.
5. Derive bracket notation, execution order, and Wigner bookkeeping from the canonical graph representation.
6. Do not store event-dependent intermediate masses in the static chain object.

## Immediate research outcomes

From `ThreeBodyDecays.jl`, the basis objects and patterns to preserve are:

- `VertexFunction` as the model for local vertex payload
- `RecouplingLS` and related recoupling types as the model for spin-coupling logic
- `Xlineshape` as the model for propagator-on-line behavior
- Wigner `d`/`D` functions as utility-level special functions with package-specific index conventions
- explicit separation between external-system data and runtime invariants

## Phase 1: Core data model

Goal: define the canonical static representation.

Deliverables:

- `AbstractVertex`
- `AbstractLineshape`
- topology container based on a line-vertex incidence matrix
- validation utilities for binary-tree cascades
- helpers for root/final/internal line classification

Questions to settle:

- exact container split between topology and payload
- explicit inclusion of the mother/root line
- storage of line quantum numbers

## Phase 2: Derived graph utilities

Goal: make the flat graph usable.

Deliverables:

- parent/child lookup tables
- topological ordering
- path-to-root queries
- bracket reconstruction utilities
- pretty-printing of topology

This phase should stay purely structural, without amplitude evaluation yet.

## Phase 3: Local physics objects

Goal: lift the prototype two-body logic into reusable package abstractions.

Deliverables:

- two-body spin helper types
- local kinematic view for one binary decay
- concrete vertex implementation modeled after `VertexFunction`
- concrete LS recoupling implementation
- first concrete propagator types or adapters

The target here is to reproduce a single local decay-step amplitude cleanly.

## Phase 4: Runtime kinematics layer

Goal: evaluate the same topology on many phase-space points without recreating the chain.

Deliverables:

- event/state object carrying external four-momenta
- derived internal line invariants
- local decay angles per vertex
- frame-transport data needed for Wigner rotations
- optional adapter to `InstructionalDecayTrees.jl` for computing these quantities from four-vectors

This phase is where the intermediate-mass issue is resolved correctly.

Recommended direction:

- keep the static `DecayChain` free of event-dependent masses and angles
- derive line invariant masses from final-state descendants
- derive vertex `(cosθ, ϕ)` values through generated instruction programs
- derive Wigner rotations through tracked instruction-path comparisons
- use helicity convention only in the first implementation
- make the `InstructionalDecayTrees.jl` bridge optional at first because it depends on unregistered four-vector tooling

## Phase 5: Cascade evaluator

Goal: automate what the current `SimpleCascade` prototype does by hand.

Deliverables:

- static overall kinematics/system object for fixed masses, spins, and line metadata
- graph-indexed amplitude input for dynamic line masses, vertex angles, and helicities
- graph-driven local-amplitude evaluation
- propagator insertion along internal lines
- summation over internal helicities
- deterministic contraction order
- optional caching of repeated local factors

The evaluator should use the topology/relation map to route local arguments:
three mass-squared values and helicities to vertices, and one invariant
mass-squared value to propagators. Static overall kinematics supplies fixed
masses and spins; the runtime input supplies dynamic invariants and angles.
Actual vertex and propagator computations should be selected by Julia dispatch.

The first supported case should be connected binary trees only.

## Phase 6: Validation against known limits

Goal: build trust early.

Deliverables:

- unit tests for topology validation
- unit tests for LS recoupling conventions
- unit tests for Wigner function conventions
- reproduction of selected `ThreeBodyDecays.jl` three-body chains as a special case
- comparison against the handwritten prototype for a small cascade

This is a key milestone before widening scope.

## Suggested first implementation order

1. Build the topology container and constructor validation.
2. Add graph utilities and bracket reconstruction.
3. Implement one concrete vertex and one concrete propagator interface.
4. Implement runtime kinematics for a small test cascade.
5. Build a minimal evaluator for one fixed example.
6. Generalize only after the example is stable.

## Near-term coding tasks

The next coding pass should likely create:

- `src/topology.jl`
- `src/lines.jl`
- `src/vertices.jl`
- `src/kinematics.jl`
- `src/evaluator.jl`
- `test/topology_tests.jl`
- `test/threebody_compat_tests.jl`

## Risks to watch

- mixing static topology with event-dependent kinematics
- over-encoding bracket structure in types instead of data
- hiding propagators inside vertices
- fixing conventions for Wigner rotations too late
- making the first evaluator too general before a minimal tree case is stable

## Recommended first success criterion

A good first milestone is:

"Represent a four-body binary cascade as a validated flat graph, reconstruct its bracket notation, derive an execution order, and evaluate one simple helicity amplitude with internal-helicity summation."

Once that works, the package architecture is on solid ground.
