# Amplitude Argument Routing

## Core correction

Do not introduce a separate "context" concept.

The problem is the split between:

- a static object for the overall kinematics/system definition
- a runtime kinematic input point
- a general amplitude call that uses topology to route local arguments

In the prototype `SimpleCascade`, every `TwoBodyDecay` already contains fixed
`TwoBodyMasses`, so local calls can read `m0`, `m1`, and `m2` from the vertex-like
object. In the general cascade this cannot work because intermediate masses are
phase-space variables.

Therefore:

- vertices contain only static information
- form factors and recouplings still require local arguments
- the `amplitude` input must provide the dynamic values needed by each local call
- the evaluator routes those input values to each vertex and propagator
- dispatch on vertex and propagator types performs the actual physics computation

## Static overall kinematics object

There should be a static object describing the overall system kinematics.

It can hold information such as:

- root/mother line id
- final-state line ids
- fixed external masses
- fixed root mass, when appropriate
- line spins and other static quantum numbers
- any static mapping between physical particles and graph line ids

This object belongs to the model and is reusable across phase-space points.

The important point is that it contains static system information, not dynamic
intermediate invariant masses.

## Static vertex object

A vertex payload should hold only objects like:

- recoupling model
- form-factor model
- static coupling labels or quantum-number constraints

It should not store:

- `m0`, `m1`, `m2`
- local angles
- helicities
- intermediate invariant masses

Those are inputs to the amplitude evaluation.

## Input shape

The runtime amplitude input should be graph-indexed.

At minimum, it should provide:

```julia
amplitude(chain, x, two_λs)
```

where:

- `x` contains runtime kinematic variables indexed by line and vertex
- `two_λs` contains helicities indexed by line

A concrete first version could be:

```julia
struct CascadeKinematics{M,A}
    line_masses2::M      # indexed by line id
    vertex_angles::A     # indexed by vertex id, values like (cosθ, ϕ)
end
```

This is not a new semantic context. It is just the kinematic input point.
It can be constructed from external four-vectors, from Dalitz variables, or from
another phase-space parameterization.

## Unified line-indexed views

The evaluator wants simple arrays:

```julia
line_masses2[line]
two_λs[line]
two_js[line]
```

but each one is assembled from mixed sources.

This should be handled before local vertex evaluation by constructing a unified
line-indexed view. The local evaluator should never branch on whether a line is
external, internal, or root.

The recommended convention is:

- line ids are canonical and already defined by the topology
- every line id has one slot in each view
- root, final, and internal lines are all included
- the source of each slot is allowed to differ

This gives uniform local code while keeping construction flexible.

## Line ordering convention

For the first implementation, prefer an explicit line-id convention:

1. final-state lines are `1:n`
2. internal lines follow
3. the root/mother line is last

This is not required by the incidence matrix in principle, but it makes input
construction and debugging much easier.

With this convention:

- external fixed data naturally fills the first `n` entries
- internal dynamic variables fill the middle entries
- root fixed or event-derived total values fill the last entry

If we later need arbitrary line ids, this convention can become a constructor
default rather than a hard requirement.

## Mass routing

For a vertex:

```text
parent -> child1 + child2
```

the evaluator obtains line ids from the topology:

```julia
l0, l1, l2 = vertex_lines(topology, vertex)
```

`incoming_line` is unambiguous because each binary vertex has exactly one
incoming line. The outgoing side needs an ordering convention: the two `+1`
entries in the incidence matrix are just children, but local physics calls need
`child1` and `child2`.

Therefore, local amplitude code should use `vertex_lines`, not raw
`outgoing_lines`. `vertex_lines` should return the outgoing children in canonical
child order, for example sorted by the smallest final-state descendant. This
makes the bracket reconstruction and the local argument routing use the same
ordering.

Then local masses are routed from the input:

```julia
m0² = x.line_masses2[l0]
m1² = x.line_masses2[l1]
m2² = x.line_masses2[l2]
```

These numbers may originate from different places:

- fixed final-state masses from the static overall kinematics object
- fixed root/mother mass from the static overall kinematics object
- internal invariant masses from the phase-space point
- external four-vectors converted by `InstructionalDecayTrees.jl`

But after the input is built, the evaluator does not care where they came from.

So `line_masses2` is best understood as a completed view:

```text
[final masses²..., internal invariant masses²..., root mass²]
```

where final/root values may be copied from static system information and internal
values come from the kinematic point.

## Form-factor call

The local form-factor call becomes:

```julia
vertex.ff(m0², m1², m2²)
```

or, if we choose mass rather than mass-squared convention:

```julia
vertex.ff(m0, m1, m2)
```

The important part is that the vertex object itself does not store these values.

Given `ThreeBodyDecays.jl` uses squared masses in form factors, the first
implementation should probably use squared masses internally:

```julia
ff(m0², m1², m2²)
```

## Recoupling call

The local recoupling call routes helicities and spins similarly:

```julia
two_λ0 = two_λs[l0]
two_λ1 = two_λs[l1]
two_λ2 = two_λs[l2]

two_j0 = line_spins[l0]
two_j1 = line_spins[l1]
two_j2 = line_spins[l2]
```

Then:

```julia
amplitude(vertex.recoupling, (two_λ1, two_λ2), (two_j0, two_j1, two_j2))
```

The spins are static line metadata from the overall system/model definition.
The helicities are amplitude input, because the evaluator will sum over internal
helicities.

So `two_js` is the simplest view:

```text
[external spins..., internal resonance spins..., root spin]
```

It belongs to the static overall kinematics/model object.

`two_λs` is more subtle because external helicities may be user-requested while
internal helicities are summation indices. The evaluator should still assemble a
complete line-indexed helicity tuple for each term in the internal-helicity sum:

```text
[external helicities..., internal running helicities..., root helicity]
```

The local vertex code then sees only `two_λs[line]`.

## Angle routing

Local angles are vertex-indexed input:

```julia
angles = x.vertex_angles[vertex]
cosθ = angles.cosθ
ϕ = angles.ϕ
```

They are computed before the amplitude call, likely by an
`InstructionalDecayTrees.jl` adapter, and then routed by vertex id.

Angles have the same mixed-source issue, but with vertex ids instead of line ids.
Each binary vertex should have exactly one local angle slot:

```julia
vertex_angles[vertex] = (cosθ = ..., ϕ = ...)
```

Those angles may be computed from external four-vectors, from a phase-space
parameterization, or from generated `InstructionalDecayTrees.jl` programs.
The evaluator should not care about the source.

The important convention is which child defines the measured direction. For the
first implementation, this should be topology-derived and deterministic, for
example:

- use the first child returned by `outgoing_lines(topology, vertex)` after
  canonical child ordering
- store the resulting angle under the vertex id

If a physics convention later needs "particle 2" directions at specific
vertices, the adapter should encode that while still returning a single
`vertex_angles[vertex]` value.

## Propagator call

Propagators are attached to internal lines. For a propagating line `l`:

```julia
σ = x.line_masses2[l]
propagator(σ)
```

The propagator stores pole parameters and model constants, while `σ` is supplied
by the input kinematic point.

## Dispatch-based computation

The general amplitude evaluator should not know the details of a specific vertex
or propagator model.

Its job is only to route local arguments:

- three mass-squared values to a vertex form factor
- local helicities and spins to a recoupling
- local angles to the Wigner/helicity factor
- one invariant mass squared to each propagator

Then Julia dispatch selects the actual computation:

```julia
amplitude(vertex::SomeVertexType, local_args...)
propagator(σ)
```

or, for split vertex payloads:

```julia
amplitude(vertex.recoupling, helicities, spins)
vertex.ff(m0², m1², m2²)
```

This means new physics models are added by defining new vertex/propagator types
and their call methods, not by changing graph traversal.

## How this fits package internals

The current structure already has most of what is needed:

- `DecayTopology` tells which line enters/leaves each vertex
- `DecayChain.vertices` stores static vertex payloads
- `DecayChain.propagators` stores static propagator payloads
- `DecayChain.propagating_lines` maps propagators to internal line ids

The missing additions are:

- a static overall kinematics/system object, especially fixed masses and `two_j`
- an input type for graph-indexed kinematic variables
- a helicity container indexed by line id
- helper functions that extract local argument tuples for a vertex
- constructors that merge static and dynamic sources into complete line-indexed
  and vertex-indexed views

## Suggested helper functions

These helpers keep routing explicit without inventing a new abstraction:

```julia
vertex_lines(topology, v) -> (l0, l1, l2)
vertex_masses2(chain, x, v) -> (m0², m1², m2²)
vertex_helicities(two_λs, topology, v) -> (two_λ0, two_λ1, two_λ2)
vertex_spins(chain, v) -> (two_j0, two_j1, two_j2)
vertex_angles(x, v) -> (cosθ, ϕ)
line_invariant(x, l) -> σ
```

Then the local amplitude implementation can be plain:

```julia
function local_amplitude(chain, x, two_λs, v)
    V = chain.vertices[v]
    two_λ0, two_λ1, two_λ2 = vertex_helicities(two_λs, chain.topology, v)
    two_j0, two_j1, two_j2 = vertex_spins(chain, v)
    m0², m1², m2² = vertex_masses2(chain, x, v)
    angles = vertex_angles(x, v)

    rec = amplitude(V.recoupling, (two_λ1, two_λ2), (two_j0, two_j1, two_j2))
    ff = V.ff(m0², m1², m2²)
    D = local_wigner_factor(two_j0, two_λ0, two_λ1 - two_λ2, angles)

    return D * rec * ff
end
```

## Static system definition

The static system can still hold known masses of external/root lines.

But these should be converted into `x.line_masses2` before or during amplitude
input construction. That way the vertex routing is uniform:

```julia
m² = x.line_masses2[line]
```

for every line, regardless of whether the value is static or dynamic.

## Important rule

The amplitude evaluator should normalize all inputs into graph-indexed arrays or
static vectors before local vertex calls.

Once inside a local vertex call, there should be no graph logic and no lookup of
global state. Just local two-body arguments.
