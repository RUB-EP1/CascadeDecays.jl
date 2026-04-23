# Open Design Questions

## Main unresolved issue: masses on intermediate lines

The key issue you pointed out is real and important:

- the object should know the particle content of the cascade
- but it should not freeze internal line masses as immutable numeric values if those masses are part of phase-space variables

## Recommended resolution

The canonical topology object should not store instantaneous intermediate masses.

Instead, distinguish clearly between:

### 1. Static line metadata

Things that belong to the line as a particle/state type:

- spin
- parity
- identity label
- which lines carry propagators
- which propagator model is attached

### 2. Runtime kinematic state

Things that vary event by event:

- four-momentum of each line
- invariant mass squared of each internal line
- local decay angles
- frame-transport information

This split is the most important design decision for the package.

## Consequence for constructors

Creating a `DecayChain` should define:

- topology
- static quantum numbers
- dynamical models

Evaluating a `DecayChain` should take:

- external event kinematics
- or a derived kinematic state for all lines

and from that compute:

- internal invariants
- local two-body systems at each vertex
- Wigner rotations and amplitudes

So the object does not need to be recreated for each phase-space point.
Only the runtime kinematic state changes.

## If a line has a pole mass

A propagating state can still carry fixed model parameters such as:

- pole mass
- nominal width
- Blatt-Weisskopf radius
- coupling constants

Those belong to the propagator model, not to the event kinematics.

So there is no contradiction:

- fixed resonance parameters live in the lineshape object
- event-dependent invariant mass is passed into that object during evaluation

## Design recommendation

Prefer this conceptual split:

- topology object: static
- model object: static
- event kinematics object: dynamic
- evaluation cache / execution plan: derived

## Questions to settle early

The next implementation pass should make decisions on:

1. Whether the root/mother line is always explicitly included in the graph.
2. Whether external lines carry their quantum numbers in a separate payload array from internal lines.
3. Whether vertex form factors are part of `AbstractVertex` or a separate payload.
4. Whether the first implementation supports only tree cascades or also non-tree diagrams.
5. Whether Wigner-rotation conventions are fixed initially to one convention or abstracted from day one.

## Strong recommendation

For the first implementation:

- support only connected binary trees
- include the root explicitly
- attach propagators to internal lines only
- treat internal invariant masses as runtime kinematics
- derive evaluation order from topology instead of storing it manually
