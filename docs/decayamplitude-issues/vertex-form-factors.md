# Make vertex form factors explicit functions of the connected masses

## Summary

`decayamplitude` currently has no first-class vertex form-factor entity. A
generic `Resonance.lineshape` callback is evaluated inside the vertex coupling
and may receive only `d1_mass` and `d2_mass`; examples typically capture a
parent mass in a Python closure. This makes the ownership, mass source, child
ordering, orbital momentum, and normalization of a form factor impossible to
inspect or serialize reliably.

A form factor belongs to a decay vertex and should have an explicit contract
based on the three lines connected by that vertex:

```text
F_v(m0, m1, m2 | form_factor_parameters)
```

where `m0` is the actual parent-line mass and `m1,m2` are the actual child-line
masses in stored vertex order. If an implementation uses squared masses, that
must be a separately declared API/schema convention rather than an implicit
choice.

Relevant current code:

- [the generic callback receives LS labels and optional daughter masses](https://github.com/KaiHabermann/decayamplitude/blob/39e0f1fb03778abb29ccab4bf09c497a2b1c34e3/src/decayamplitude/resonance.py#L147-L187)
- [`DecayChainNode` derives the two child masses at the connected vertex](https://github.com/KaiHabermann/decayamplitude/blob/39e0f1fb03778abb29ccab4bf09c497a2b1c34e3/src/decayamplitude/chain.py#L167-L212)

## Proposed ownership model

Keep the factors separate:

```text
internal line: X_R(sigma_R | propagator parameters)
vertex:        H_v * F_v(m0,m1,m2 | FF parameters) * D_v*
```

The vertex owns:

- the ordered parent, child 1, and child 2 line references;
- its helicity or LS coupling terms;
- the orbital momentum relevant to each form-factor term;
- the form-factor analytic type and radius/scale parameters;
- an explicit normalization policy.

The propagator may own nominal decay-channel masses and channel `l,s` metadata
for a running-width model, but it must not receive the vertex's three event
masses or active LS term at evaluation time.

## Normalization must be data

For a Blatt-Weisskopf-type factor, these are different models:

```text
raw:              F_L(q)
pole-normalized:  F_L(q) / F_L(q0)
```

The API and serialized model must therefore record a normalization enum such
as `raw`, `pole`, or `explicit_reference`, together with the definition of
`q0`. It must also distinguish ordinary integer `L` from doubled-spin `2L`.

## Mass policies

The following quantities should not be conflated:

- event parent mass `m0` and child masses `m1,m2`, used by a dynamic vertex
  factor;
- nominal particle/channel masses, used for threshold and pole-reference
  definitions;
- the pole mass in an internal-line propagator denominator.

In particular, an external particle's reconstructed four-vector mass must not
silently replace the pole mass in a Breit-Wigner denominator.

## Suggested acceptance criteria

- [ ] A vertex form factor is a first-class, inspectable object separate from a
      propagator.
- [ ] Its evaluation interface receives exactly the three connected masses (or
      a named `VertexMasses` value) plus form-factor parameters.
- [ ] Parent, child 1, and child 2 ordering is preserved and tested.
- [ ] Mass-valued versus squared-mass-valued inputs are explicit.
- [ ] Ordinary `L`, scale/radius units, momentum definition, and normalization
      policy are serializable.
- [ ] Raw and pole-normalized Blatt-Weisskopf models have distinct tests.
- [ ] Nested-cascade tests verify that off-shell internal child masses are
      routed to the correct vertex.
- [ ] Propagator callbacks no longer receive vertex masses or vertex LS labels.

This issue complements the separate proposal to make propagators independent
internal-line entities evaluated only as `X(sigma | parameters)`.
