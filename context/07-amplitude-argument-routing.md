# Amplitude Argument Routing

Topology indices (`line`, vertex index, `address`, incidence matrix) are documented in
[`docs/src/notation.md`](../docs/src/notation.md). This note covers only how runtime
inputs are routed into local amplitude calls.

## Naming

| Name | Meaning |
|------|---------|
| `line` | Topology line id (`1:nlines`) |
| `vertex_ind` | Topology vertex index (`1:nvertices`) |
| `address` | Bracket key in `address => payload` constructor pairs |
| `vertex_func` | Static `VertexFunction` at one vertex: `chain.vertices[vertex_ind]` |
| `propagator` | Static lineshape callable: `chain.propagators[i]` on `propagating_lines[i]` |

Do not use bare `vertex` for the func — the word also means the graph slot, the index,
and the `vertices = (...)` keyword of address/func pairs.

## Core correction

Do not introduce a separate "context" concept.

The problem is the split between:

- a static object for the overall kinematics/system definition
- a runtime kinematic input point
- a general amplitude call that uses topology to route local arguments

In the prototype `SimpleCascade`, every `TwoBodyDecay` already contains fixed
`TwoBodyMasses`, so local calls can read masses from the decay object. In the general
cascade, intermediate masses are phase-space variables.

Therefore:

- `vertex_func` holds only static information (recoupling, form factor)
- form factors and recouplings still require local arguments at evaluation time
- the amplitude input provides dynamic values for each local call
- the evaluator routes those values to each `vertex_func` and `propagator`
- dispatch on `VertexFunction` and propagator types performs the physics

## Static system object

`CascadeSystem` holds reusable static information: external and root masses/spins,
line-to-particle mapping, and related quantum numbers. It does not store dynamic
internal invariant masses.

## Static vertex func

A `vertex_func` (`VertexFunction`) holds:

- recoupling model (`vertex_func.h`)
- form-factor model (`vertex_func.ff`)

It should not store `m0`, `m1`, `m2`, local angles, helicities, or intermediate masses.
Those are routed in at evaluation time.

Constructor input is `address => vertex_func`; storage is `chain.vertices[vertex_ind]`.

## Input shape

```julia
amplitude(chain, system, x)
amplitude(chain, system, x, external_two_λs)
```

- `x::CascadeKinematics` — runtime kinematics indexed by `line` and vertex index
- `external_two_λs` — optional external helicity selection

```julia
struct CascadeKinematics{...}
    line_masses2::SVector{Nl,T}   # indexed by line
    vertex_angles::SVector{Nv,A}  # indexed by vertex_ind
end
```

This is the kinematic input point, not a new semantic context. It can be built from
four-vectors (`cascade_kinematics`), Dalitz variables, or other parameterizations.

## Unified line-indexed views

Before local evaluation, assemble complete line-indexed views:

```julia
line_masses2[line]
two_λs[line]
two_js[line]   # from line_two_js(chain, system)
```

Every line id has one slot. The evaluator does not branch on final vs internal vs root
when reading these arrays.

Construction may mix static system data and dynamic `x` values; after assembly, local
code only indexes by `line`.

## Mass routing

For vertex index `vertex_ind`:

```julia
l0, l1, l2 = vertex_lines(chain, vertex_ind)
m0², m1², m2² = vertex_masses2(chain, x, vertex_ind)
# equivalent to x.line_masses2[l0], x.line_masses2[l1], x.line_masses2[l2]
```

Use `vertex_lines`, not raw `outgoing_lines`, so child order matches bracket
addresses and physics conventions.

## Form-factor and recoupling calls

```julia
vertex_func.ff(m0², m1², m2²)
ThreeBodyDecays.amplitude(vertex_func.h, (two_λ1, two_λ2), (two_j0, two_j1, two_j2))
```

Squared masses match `ThreeBodyDecays.jl` form-factor conventions.

Helicities come from the line-indexed `two_λs`; spins from `line_two_js(chain, system)`.
Internal line helicities are summed later; external ones may be fixed by the caller.

## Angle routing

```julia
angles = vertex_angles(x, vertex_ind)
# angles.cosθ, angles.ϕ
```

One angle slot per vertex index, usually from generated `InstructionalDecayTrees.jl`
programs (`cascade_kinematics`, `helicity_angle_programs`).

## Propagator call

For internal line `line` with propagator slot `i`:

```julia
σ² = line_invariant(x, line)
chain.propagators[i](σ²)
```

## Dispatch-based computation

The graph traversal only routes arguments. Physics is in typed callables:

```julia
routed_vertex_amplitude(vertex_func::VertexFunction, masses2, helicities, spins, angles)
propagator(σ²)
```

New models extend `VertexFunction` / propagator types and their methods, not traversal.

## Implemented routing helpers

```julia
vertex_lines(topology, vertex_ind) -> (l0, l1, l2)
vertex_masses2(chain, x, vertex_ind) -> (m0², m1², m2²)
vertex_helicities(chain, two_λs, vertex_ind) -> (two_λ0, two_λ1, two_λ2)
vertex_spins(chain, system, vertex_ind) -> (two_j0, two_j1, two_j2)
vertex_angles(x, vertex_ind) -> (cosθ, ϕ)
line_invariant(x, line) -> σ²
```

Local evaluation pattern:

```julia
function routed_vertex_amplitude(vertex_func, masses2, helicities, spins, angles)
    # Wigner rotation × recoupling(vertex_func.h, ...) × vertex_func.ff(masses2...)
end
```

## Important rule

Normalize inputs into graph-indexed arrays before local calls. Inside
`routed_vertex_amplitude`, there is no graph logic — only local two-body arguments.
