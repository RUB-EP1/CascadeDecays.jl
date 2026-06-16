# Lessons learned from the three-body Wigner-rotation debugging

This debugging session exposed several real bugs, but the more important lesson is
architectural: the current API makes it too easy to combine a model, a topology,
a reference plane, and kinematic variables that are not mutually consistent.
When that happens, the code still produces numbers, but the numbers no longer
represent a well-defined amplitude convention.

## What went wrong

### Reference topology was implicit where it had to be explicit

For a three-body amplitude, the reference topology used to define the production
orientation must match the reference topology used by the amplitude evaluation.
We had places where `ref_k` was passed separately from the model construction
and from the kinematic-point construction. That made it possible to evaluate a
chain in one convention and compare it against angles computed in another.

This was not a numerical issue. It was a representational issue: the relevant
convention lived in several disconnected arguments, so the type system and the
call structure could not protect us.

### Ordered topologies were accidentally canonicalized

We discovered that `DecayTopology(((3,1),2))` was being canonicalized to
`DecayTopology(((1,3),2))`. That is a bug because the child order is physical.
The two descriptions have the same invariant mass pair, but they do not define
the same helicity convention, local axes, Wigner paths, or spin phases.

This is a good example of a mathematical simplification that is valid for some
invariant quantities but invalid for amplitudes. A topology object for amplitude
work must preserve ordered structure, not merely the unordered decay tree.

### Kinematics and amplitudes could be built with incompatible conventions

`KinematicPoint` and the cascade amplitude machinery were allowed to choose
their topology/reference information independently. That made it possible to ask
a syntactically valid but physically inconsistent question:

```julia
angles = KinematicPoint(..., topology_a)
amplitude(..., reference_topology_b, angles)
```

The library should make the consistent path the easiest path. Ideally, the model
or evaluation context should carry the reference topology, child ordering, and
frame convention together, and downstream computations should derive from that
single source of truth.

### Frame handling hid a rest-frame edge case

The initial-frame choice matters. We found that applying a helicity-root-frame
boost to a system that is already numerically in its rest frame can introduce
spurious behavior. This is especially dangerous because the event may be
conceptually in the rest frame while floating-point noise makes the code follow a
different branch.

The fix was to guard this case and fall back to the current frame when the parent
is already effectively at rest. The broader lesson is that frame conventions
should be explicit, tested, and attached to the kinematic construction, not left
as a hidden procedural detail.

### Wigner rotations need branch-aware tests

The `Λ + K* + Δ` comparison showed that the ThreeBodyDecays convention for the
Δ contribution corresponds to `Λ + K* - Δ` in the CascadeDecays convention. This
is not just an arbitrary sign: it follows from an over-`2π` Wigner-rotation
branch, which flips spin-1/2 states.

This is exactly the kind of issue that can be invisible in isolated chain tests.
A single chain may match up to a global sign, but coherent sums are sensitive to
relative phases. Therefore, tests must include coherent multi-chain sums, not
only per-chain absolute checks.

### Singular or identical Wigner paths need explicit semantics

When the reference path and target path are identical, the Wigner rotation should
be exactly trivial. Letting the generic Euler-angle decomposition handle that
case exposed branch and singularity sensitivity. Treating identical paths as a
special case is not a hack; it encodes the intended geometry directly.

## What this says about the architecture

The code currently separates concepts that should probably travel together:

- the decay topology;
- the ordered child convention;
- the reference topology or reference particle;
- the production-plane orientation;
- the kinematic-point topology;
- the initial-frame convention;
- the Wigner-rotation path.

Because these are separate knobs, users can accidentally build a physically
incoherent evaluation. The library then faithfully evaluates the incoherent
request instead of rejecting it or making it difficult to express.

The architecture should move toward an explicit evaluation context, or an
equivalent object, that binds these choices together. A model should know the
topology and reference convention it is written in. Kinematic variables should be
computed from that same convention unless the user intentionally requests a
conversion. Amplitude evaluation should receive one coherent context rather than
a loose collection of compatible-looking arguments.

## Concrete guardrails to add

1. Preserve ordered topology everywhere.
   Never canonicalize child order in amplitude-facing topology constructors.
   If an unordered invariant-pair object is useful, it should be a separate type.

2. Make reference conventions first-class.
   `ref_k`, reference topology, and production orientation should be stored in a
   shared model or evaluation context rather than passed independently at several
   call sites.

3. Validate compatibility before evaluation.
   If a `KinematicPoint` was built for topology `(12)3`, evaluating it as if it
   belonged to `(23)1` or `(31)2` should require an explicit conversion or fail
   loudly.

4. Separate invariant topology from helicity topology.
   The invariant pair `(i,j)` determines masses and angles such as
   `cosθ_ij`, but the ordered pair `(i,j)` determines axes and spin phases. The
   API should reflect this distinction.

5. Add coherent-sum regression tests.
   Tests should cover `Λ`, `K*`, and `Δ` individually, but also coherent sums
   such as `Λ + K*`, `Λ + K* + Δ`, and the corresponding ThreeBodyDecays
   convention `Λ + K* - Δ`.

6. Test both current-frame and boosted-frame construction.
   The same physical event should give the same kinematic point when constructed
   in `CurrentFrame()` and in a properly handled helicity-root frame, including
   the numerically-rest-frame edge case.

7. Make branch conventions visible in documentation.
   Spin phases from Wigner rotations are convention-dependent. The documentation
   should state when two systems agree directly and when they agree only after a
   known phase convention, such as the Δ-chain sign in the TBD comparison.

## A design direction

The ideal user-facing workflow should look less like assembling independent
pieces and more like asking a model for its own kinematics:

```julia
ctx = EvaluationContext(model; reference = DecayTopology(((1,2),3)),
                        frame = CurrentFrame())
point = kinematic_point(ctx, event)
amp = amplitude(ctx, point)
```

The exact names do not matter. The important part is ownership: the context owns
the topology, child order, reference plane, and frame convention. That way the
same object that computes the kinematic variables is also the object that knows
how the amplitude will interpret them.

This would not eliminate all Wigner-rotation subtlety, but it would move many
mistakes from "valid code with wrong physics" to either "impossible to write" or
"explicit conversion requested here." That is the architecture we should aim for.
