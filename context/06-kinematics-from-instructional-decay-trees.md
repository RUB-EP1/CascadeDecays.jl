# Kinematics From `InstructionalDecayTrees.jl`

This note describes the planned optional kinematics adapter. The current package
implementation only contains the static topology/system objects and the
graph-indexed routing scaffold.

## Why this package fits

`InstructionalDecayTrees.jl` is a good match for the runtime-kinematics layer of
`CascadeDecays.jl`.

It provides a declarative DSL for:

- boosting to helicity frames
- using the particle-2 helicity convention
- measuring invariant masses
- measuring `(cosθ, ϕ)` or `(m, cosθ, ϕ)`
- tracking accumulated Lorentz transformations
- extracting relative Wigner ZYZ angles

This is exactly the part we should not hard-code into the static cascade model.

## Responsibility split

`CascadeDecays.jl` should own:

- flat topology
- line and vertex payloads
- amplitude factors
- helicity sums
- contraction order

`InstructionalDecayTrees.jl` should be used to provide:

- event-dependent internal invariant masses
- local vertex decay angles
- frame paths for each line or vertex
- relative Wigner rotations between helicity-frame paths

So the integration should be an adapter layer, not a replacement for our topology container.

## Static versus runtime data

The `DecayChain` object remains static:

- graph topology
- propagator models
- vertex models
- static quantum numbers

The event kinematics object is dynamic:

- external four-vectors
- internal line four-vectors
- internal invariant masses
- local decay angles
- Wigner rotations

This resolves the intermediate-mass issue: intermediate masses are computed per event from four-vectors and are not stored in the static chain.

## Proposed object model

Introduce a runtime object similar to:

```julia
struct EventKinematics{N, M, A, W}
    external::N
    line_masses2::M
    vertex_angles::A
    wigner_angles::W
end
```

where:

- `external` is the tuple of input four-vectors
- `line_masses2[line]` gives invariant mass squared for each graph line
- `vertex_angles[vertex]` gives `(cosθ, ϕ)` for the selected child direction in the local rest frame
- `wigner_angles` stores relative ZYZ rotations needed to compare helicity frames

The exact container can be static arrays, named tuples, or a hybrid. The first implementation should prioritize correctness and clear semantics.

## How to derive masses

For each graph line, compute its final-state descendants from the topology.

Then for a line `ℓ`:

```julia
indices = final_descendants(topology, ℓ)
MeasureInvariant(tag_for_line(ℓ), indices)
```

or compute directly from external four-vectors if the backend exposes the summed four-vector helper.

This gives:

- final-state line mass squared
- internal line invariant mass squared
- root invariant mass squared

These values feed:

- propagators on internal lines
- local two-body systems at vertices
- form factors that depend on parent and daughter masses

## How to derive local angles

For each binary vertex:

```text
parent -> child_a + child_b
```

the local angle program is:

1. go to the parent rest frame
2. measure one child direction as `(cosθ, ϕ)`
3. optionally go deeper following the selected branch

In `InstructionalDecayTrees.jl`, this is expressed with:

```julia
(
    ToHelicityFrame(parent_final_indices),
    MeasureCosThetaPhi(tag_for_vertex(v), child_final_indices),
)
```

For the initial implementation, the adapter should use helicity convention only.
Where the Jacob-Wick particle-2 step is required by the helicity path, the generated
program should use `ToHelicityFrameParticle2`; this is part of the fixed helicity
construction rather than a user-facing convention switch.

## How to derive Wigner rotations

`InstructionalDecayTrees.jl` provides:

- `TrackedState`
- `compare_instruction_paths`
- `wigner_zyz`
- `wigner_zyz_so3`
- `wigner_zyz_su2`

For `CascadeDecays.jl`, the important derived quantity is the relative frame rotation between:

- the natural path of a local decay chain
- the selected global/reference helicity basis

The adapter should generate instruction paths from graph paths and compare them with:

```julia
cmp = compare_instruction_paths(path_reference, path_other, external_fourvectors)
w = wigner_zyz(cmp.relative)
```

Then `w` provides `(ϕ, θ, ψ)` for Wigner `D`/`d` factors.

## Suggested API shape

An initial public API could look like:

```julia
kinematics(chain, external_fourvectors)
```

returning an `EventKinematics`.

Internally it would:

1. compute final descendants for every graph line
2. build `InstructionalDecayTrees.jl` measurement programs
3. execute the programs on the input four-vectors
4. store masses and angles in graph-indexed containers

The amplitude evaluator should then accept:

```julia
amplitude(chain, kin, helicities)
```

instead of accepting raw masses and angles directly.

## Dependency strategy

Because `InstructionalDecayTrees.jl` depends on an unregistered `FourVectors.jl`, it may be better to make it an optional integration first:

- keep core topology/evaluator independent
- add an extension or adapter file later
- provide a documented interface that any kinematics backend can implement

This keeps `CascadeDecays.jl` usable even before all external dependencies are registered.

## First implementation target

The first kinematics milestone should be:

1. Given a `DecayTopology` and external four-vectors, compute `line_masses2`.
2. For the topology `((1,2),3)`, compute local vertex angles for:
   - root vertex: direction of `(1,2)` or `3` in the root frame
   - isobar vertex: direction of `1` in the `(1,2)` frame
3. Verify against a direct `InstructionalDecayTrees.jl` program.

Only after this should we wire these values into amplitude evaluation.
