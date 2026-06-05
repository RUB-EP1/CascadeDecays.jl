# Internal notation {#notation}

Short reference for indices and the flat topology used inside `CascadeDecays.jl`.

## Names

| Name | Range | Meaning |
|------|-------|---------|
| `line` | `1:nlines` | Row in the incidence matrix: a particle line (final, internal, or root) |
| vertex index | `1:nvertices` | Column in the incidence matrix: a binary decay vertex |
| `address` | — | Bracket key such as `(2, 3)`; resolves via [`vertex_for`](@ref) / [`line_for`](@ref) |
| `vertex_func` | — | Local amplitude model (`VertexFunction`): `chain.vertices[v]` for vertex index `v` |
| `propagator` | — | Lineshape payload, **not** a line id: entry in `chain.propagators` |

User input uses `address => vertex_func` pairs in the `vertices` keyword; [`DecayChain`](@ref) stores the
funcs in [`vertex_address`](@ref) / root-first preorder as `chain.vertices`.

Counts for a tree with `n` final-state particles:

- `nfinal = n`
- `nvertices = n − 1`
- `nlines = nfinal + nvertices` (finals + internal lines + root)
- `Np = nvertices` (one propagator per internal line)

Bracket addresses such as `(2, 3)` or `((1, (2, 3)), 4)` are user-facing keys. They resolve to a `line` or vertex index via [`line_for`](@ref) and [`vertex_for`](@ref).

## Example

```math
B^{+} \to \psi\, K^{+}\, K^{-}\, \pi^{+}
```

Decay chain:

```text
B+  →  X0  K+
X0  →  ψ   K*0
K*0 →  K−  π+
```

Final-state labels `1…4` and bracket tree (right child = spectator at each step):

```text
        B+  (root)
       /  \
     X0    K+        finals: 1=ψ, 2=K−, 3=π+, 4=K+
    /  \
   ψ   K*0
      /  \
     K−   π+
```

```julia
topology = DecayTopology(((1, (2, 3)), 4))
# bracket(topology) == "((1,(2,3)),4)"
```

### Line numbering

| `line` | Role | Particle / system |
|-------:|------|-------------------|
| 1–4 | final | ψ, K−, π+, K+ |
| 5 | internal | K*0 (from `K−` + `π+`) |
| 6 | internal | X0 (from ψ + K*0) |
| 7 | root | B+ |

Internal lines are assigned in **postorder**; the root line is always `nlines`.

### Vertex indices

Vertex indices are the **columns** of the incidence matrix, in **root-first preorder** (same order as `DecayChain` `vertices`):

| vertex index | Bracket address | Decay |
|-------------:|-----------------|-------|
| 1 | `((1, (2, 3)), 4)` | B+ → X0 + K+ |
| 2 | `(1, (2, 3))` | X0 → ψ + K*0 |
| 3 | `(2, 3)` | K*0 → K− + π+ |

At each vertex index, [`vertex_lines`](@ref) returns `(parent, child1, child2)` as line ids. For example index 2: parent line 6 (X0), children lines 1 (ψ) and 5 (K*0).

### Incidence matrix

`relation[line, v]` with vertex index `v` uses the sign convention documented on [`DecayTopology`](@ref):

- `−1` — line enters the vertex (incoming parent)
- `+1` — line leaves the vertex (outgoing child)
- `0` — no incidence

Rows = lines, columns = vertex indices:

```text
          v1  v2  v3
  L1 ψ    0   1   0
  L2 K−   0   0   1
  L3 π+   0   0   1
  L4 K+   1   0   0
  L5 K*0  0   1  −1
  L6 X0   1  −1   0
  L7 B+  −1   0   0
```

### Propagators vs lines

Only **internal** lines (5 and 6 here) carry propagators. A propagator payload sits in `chain.propagators[i]`; the corresponding line is `propagating_lines(chain)[i]`:

| `i` | `line` | Resonance |
|----:|-------:|-----------|
| 1 | 5 | K*0 |
| 2 | 6 | X0 |

Final lines (1–4) and the root (7) are lines without an attached propagator.

### Amplitude axes

[`line_amplitude_tensor`](@ref) builds a buffer with **one axis per line**, sized by the spin/helicity multiplicity on that line. [`external_helicity_amplitude`](@ref) sums over internal line axes ([`propagating_lines`](@ref)); the result has one axis per external line (finals in [`finallines`](@ref) order, plus the root helicity last).
