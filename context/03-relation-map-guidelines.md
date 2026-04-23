# Relation Map Guidelines

## Canonical meaning

The `relation` map should be treated as a line-vertex incidence matrix.

Rows are lines.
Columns are vertices.

Each vertex is a binary decay coupling.

Each line is a particle/state flowing between vertices.

## Which lines are included

The matrix should include:

- all final-state particles
- all intermediate particles
- optionally the mother/root line

Including all lines is strongly preferred because it makes the representation self-contained.

## Sign convention

A simple and useful convention is:

- `-1`: line enters the vertex
- `+1`: line leaves the vertex
- `0`: line not attached to the vertex

For a binary decay tree, each column then contains exactly:

- one `-1`
- two `+1`

If a different convention is later chosen, it should preserve this same information content.

## Tree constraints

For ordinary cascade decays without loops:

- every non-root internal line is produced exactly once and consumed exactly once
- every final-state line is produced exactly once and never consumed
- the root is consumed once and never produced

These rules allow fast validation of the topology object.

## Example

For the bracket structure `((1,2),3)`, including the mother line `0` and intermediate line `R12`, one possible line ordering is:

- line 1: particle 1
- line 2: particle 2
- line 3: particle 3
- line 4: `R12`
- line 5: mother `0`

and two vertices:

- `v1`: `0 -> R12 + 3`
- `v2`: `R12 -> 1 + 2`

Then the incidence matrix can be:

```text
        v1  v2
1        0  +1
2        0  +1
3       +1   0
4       +1  -1
5       -1   0
```

This flat data is enough to derive:

- the bracket topology
- parent-child relations
- leaf order
- evaluation order

## Bracket reconstruction

The matrix should be the canonical topology representation.
Bracket notation should be derived from it, not stored independently.

That avoids duplicated topology state.

## Propagator attachment rule

Propagators should be attached only to internal lines.

Not every line has a propagator:

- final-state lines do not
- the root line usually does not
- only intermediate resonant or exchanged lines do

So the package needs a mapping:

- graph line id
- internal propagator payload index

This should be explicit and not inferred from row position alone.

## Computational ordering

Evaluation order should also be derived, not user-entered.

For a tree decay this usually means:

1. order vertices by depth from leaves to root, or vice versa depending on the chosen contraction convention
2. infer which propagator factors appear between which local amplitudes
3. derive frame-transport paths from ancestor relations

This is one of the main benefits of the flat graph approach.

## Validation checklist for constructor

The topology constructor should validate at least:

- every vertex has one incoming and two outgoing lines
- there is exactly one root if roots are stored
- final-state lines are leaves
- no disconnected components
- no cycles
- the number of vertices is compatible with a binary tree

For a connected binary tree with root included:

- `Nv = Ninternal + 1`
- `Nlines = Nfinal + Ninternal + 1`

Those identities are useful sanity checks.
