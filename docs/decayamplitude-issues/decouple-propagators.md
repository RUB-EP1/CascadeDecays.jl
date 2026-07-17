# Decouple propagators from vertex LS/helicity couplings

## Summary

The propagating resonance and its production/decay vertices are currently
represented by the same `Resonance` object. In particular,
`Resonance.helicity_from_ls` evaluates `self.lineshape` inside the LS-coupling
sum and calls it with

```python
self.lineshape(two_l, two_s, *params, d1_mass=..., d2_mass=...)
```

while the direct-helicity path calls the same callback with
`(two_h1, two_h2)`.

This couples an internal-line propagator to the particular vertex term used to
produce or decay it. A propagator should instead be its own entity, evaluated
from the actual invariant mass only:

```text
X_R(sigma_R | propagator_parameters)
```

Here `sigma_R = p_R^2` (or an explicitly documented mass-valued equivalent).
No vertex helicities, LS labels, or event-level daughter masses should be
passed at evaluation time.

Relevant current code:

- [LS path in `resonance.py`](https://github.com/KaiHabermann/decayamplitude/blob/39e0f1fb03778abb29ccab4bf09c497a2b1c34e3/src/decayamplitude/resonance.py#L147-L187)
- [direct-helicity path in `resonance.py`](https://github.com/KaiHabermann/decayamplitude/blob/39e0f1fb03778abb29ccab4bf09c497a2b1c34e3/src/decayamplitude/resonance.py#L233-L249)
- [child masses supplied by `DecayChainNode`](https://github.com/KaiHabermann/decayamplitude/blob/39e0f1fb03778abb29ccab4bf09c497a2b1c34e3/src/decayamplitude/chain.py#L167-L212)

## Proposed ownership model

Introduce a first-class propagator/internal-line object. It may own a
declarative list of decay-channel definitions, with each channel carrying its
own quantum numbers and parameters independently, for example:

```python
DecayChannel(
    two_l=...,
    two_s=...,
    particle1=...,
    particle2=...,
    parameters=...,
)
```

Those channel definitions belong to the propagator model; they must not be
inferred from, or injected by, the currently evaluated vertex coupling.
Evaluation remains a function only of the line invariant and propagator
parameters:

```python
propagator(sigma, parameters)
```

The vertices separately supply helicity/LS couplings, angular functions, and
vertex form factors. Schematically,

```text
amplitude = product_R X_R(sigma_R | p_R)
            * product_v V_v(connected momenta, vertex parameters)
```

This also makes serialization unambiguous: internal lines, their channel
tables, and their parameters can be serialized separately from vertices and
their couplings.

## Why this matters

1. Multiple vertex LS terms that share one internal resonance should share one
   propagator value, rather than re-evaluating a nominal "line shape" inside
   every term.
2. The propagator cannot accidentally depend on particle ordering, helicity
   convention, or vertex-specific LS recoupling.
3. Event-dependent daughter masses cannot leak into a Breit-Wigner denominator
   through a generic callback. The denominator depends on the actual
   propagating invariant; nominal channel masses belong to channel/model
   parameters.
4. A root decay has a vertex but normally no parent propagator, which becomes
   explicit in the model.
5. The split provides a stable serialization boundary compatible with other
   cascade-amplitude frameworks.

## Suggested acceptance criteria

- [ ] A propagator/internal-line is represented independently of a decay vertex.
- [ ] Its public evaluation interface accepts only the actual line invariant
      and propagator parameters.
- [ ] Decay-channel metadata, including separate `l` and `s`, is stored on the
      propagator/channel model and is not taken from the active vertex term.
- [ ] A chain computes and evaluates each internal propagator once per
      event/line.
- [ ] Vertex helicities, LS labels, and event-level daughter masses are never
      passed to the propagator callback.
- [ ] The root node does not acquire a propagator merely because it has a decay
      vertex.
- [ ] Tests cover one resonance shared by multiple LS vertex terms and show a
      common propagator factor.
- [ ] A migration path is documented for the current
      `lineshape(l, s, ..., d1_mass, d2_mass)` callback API.

This issue is about the entity boundary and evaluation contract. Vertex
form-factor ownership is a separate, complementary concern.
