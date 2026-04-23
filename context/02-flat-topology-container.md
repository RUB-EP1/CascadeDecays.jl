# Flat Topology Container

## Target idea

The intended canonical container is a flat graph:

```julia
struct DecayChain{n, Np, Nv, P<:AbstractLineshape, V<:AbstractVertex}
    Pv::SVector{Np, P}
    Vs::SVector{Nv, V}
    relation::SMatrix{Np+n, Nv, Int}
end
```

with:

- `n`: number of final-state particles
- `Np`: number of propagating internal lines
- `Nv`: number of binary decay vertices

This direction is correct in spirit, but the meaning of the dimensions needs one refinement.

## Recommended interpretation

The relation matrix should cover all lines, not only propagating ones.

So the conceptual object is:

- rows: lines
- columns: vertices

where lines include:

- final-state external lines `1, ..., n`
- internal propagating lines
- optionally the root / mother line

Then each column represents one binary decay:

- one incoming line
- two outgoing lines

This is a line-vertex incidence matrix, not a propagator-vertex map.

## Why this is better than recursive nesting

This keeps the package aligned with Julia strengths:

- concrete flat storage
- static array friendliness
- easier type stability
- easier serialization
- easier graph traversals

It also avoids a recursive type tree like:

`Decay(Decay(Decay(1,2),3),4)`

which is elegant for notation but awkward as the canonical computational representation.

## Recommended structural split

The container should probably be split conceptually into three layers.

### 1. Topology

Pure connectivity:

- line ids
- vertex ids
- incidence map
- root line
- final-state line ids

### 2. Payload

Physics attached to graph elements:

- propagators attached to internal lines only
- vertices attached to vertices
- quantum numbers attached to lines and possibly vertices

### 3. Derived execution data

Computed from topology plus conventions:

- evaluation order
- parent/child lookup tables
- line depths
- path-to-root tables
- bracket reconstruction
- frame transport / Wigner-rotation bookkeeping

This separation will make the package easier to evolve.

## Recommended revision of the prototype shape

The current prototype puts only propagators and vertices in typed arrays and lets `relation` carry the rest.
That is close, but one more explicit line-level payload layer would help.

For example:

```julia
struct DecayTopology{Nl,Nv}
    relation::SMatrix{Nl,Nv,Int}
    root::Int
    finals::SVector
end

struct DecayDynamics{Np,Nv,P<:AbstractLineshape,V<:AbstractVertex}
    propagators::SVector{Np,P}
    vertices::SVector{Nv,V}
    propagating_lines::SVector{Np,Int}
end
```

Then a higher-level chain object can combine them.

This is not the only valid split, but it keeps the graph semantics explicit.

## Important invariant

For a binary decay cascade, each vertex column should represent:

- exactly one incoming line
- exactly two outgoing lines

If the root line is included, then:

- the root appears as incoming only once
- each final-state line appears as outgoing only once and never incoming again
- each internal line appears as outgoing once and incoming once

This gives a clean tree structure while still being stored as a flat matrix.
