# Four-body integration: `B⁺ → ψ(4040) K⁺`, `ψ(4040) → Dₓ⁺ D⁻`, `Dₓ⁺ → D⁰ π⁺`

This literate example demonstrates how the current `CascadeDecays.jl` internals
route kinematic variables through a flat cascade topology.

The input event is taken from:

- <https://github.com/RUB-EP1/B2DxDK.jl/blob/fe8f63cad0a64c5d757cd01d7f61be0644026b28/data/crosscheck_event.json>

The angle program follows the convention used in:

- <https://github.com/RUB-EP1/B2DxDK.jl/blob/fe8f63cad0a64c5d757cd01d7f61be0644026b28/notebooks/all_angles.jl>

We use helicity convention only. The goal here is not yet a production physics
model, but a clean integration point: static graph and model payloads live in
`CascadeDecays.jl`, while event-dependent masses and angles are computed from
four-vectors and routed into local vertex and propagator calls.

````julia
using CascadeDecays
using FourVectors
using HadronicLineshapes
using StaticArrays
using ThreeBodyDecays: VertexFunction, RecouplingLS

format_instruction(inst) = string(nameof(typeof(inst)), "(", join((repr(getfield(inst, i)) for i in 1:fieldcount(typeof(inst))), ", "), ")")
````

````
format_instruction (generic function with 1 method)
````

## 1. Event four-vectors

The final-state particle labels are chosen to match the
`InstructionalDecayTrees.jl` angle notebook:

```text
1: D⁰
2: π⁺
3: D⁻
4: K⁺
```

````julia
pDminus = FourVector(0.8634762475578601, -0.2273640501540901, 0.5897962254778486; E = 2.1542368373711818)
pD0 = FourVector(-0.41098561980779524, 0.4602903155362023, 0.1331926788204417; E = 1.9687832721903593)
pKplus = FourVector(-0.4483316412801869, -0.23415599517164887, -0.7305348220308255; E = 1.0164784292725653)
piplus = FourVector(-0.004158986333048991, 0.0012297296374225927, 0.00754591768614862; E = 0.13984149613322175)

objs = (pD0, piplus, pDminus, pKplus);

P_Dx = pD0 + piplus
P_psi = P_Dx + pDminus
P_B = P_psi + pKplus;
````

The event record contains a few cross-check values. We keep them close to the
source values to guard against accidental reordering.

````julia
reference_mψ2 = 17.3825
reference_cosθ_D_in_Dx = 0.8746538492596707
reference_ϕ_D_in_Dx = -2.84901364039537;

@assert isapprox(mass(P_psi)^2, reference_mψ2; atol = 5e-5)
````

## 2. Topology from bracket notation

The bracket notation `(((1,2),3),4)` defines the full binary cascade:

```text
(((D⁰,π⁺),D⁻),K⁺)
```

`DecayTopology` converts this notation into the flat line-vertex incidence
matrix. The generated convention is final-state lines first, internal lines
next, and root last.

````julia
topology = DecayTopology((((1, 2), 3), 4))

@assert has_canonical_line_order(topology)
@assert bracket(topology) == "(((1,2),3),4)"
````

## 3. Static system information

The static system stores only external/root spins and fixed external/root
masses. Internal spins belong to the propagator specs below, and internal
masses are supplied by the event kinematics.

We use pseudoscalar external particles and a pseudoscalar `B⁺`. Spins are
represented as doubled integers.

````julia
two_js = (
    0, # D⁰
    0, # π⁺
    0, # D⁻
    0, # K⁺
    0, # B⁺
)

masses2 = (mass.(objs) .^ 2..., mass(P_B)^2)
system = CascadeSystem(two_js, masses2);
````

## 4. Runtime kinematic input from the topology

`CascadeKinematics` is the completed input point from which the package can
route local invariants and angles by topology address.

The angle programs are generated directly from the topology. It is useful to
look at them explicitly: each entry is the sequence of frame changes and one
final `MeasureCosThetaPhi(...)` instruction for one binary decay in the
cascade.

In helicity convention:

- `ToHelicityFrame((...))` boosts to the rest frame of the addressed subsystem
- the final `MeasureCosThetaPhi(...)` returns the local helicity angles
  `(cosθ, ϕ)` for the first child of that binary decay

The initial configuration of four-vectors is assumed to be already aligned so
that the subsequent helicity-angle measurements are performed in the desired
convention.

After that, `cascade_kinematics(...)` evaluates the generated programs on the
event four-vectors and assembles the full kinematic input object.

Here we inspect the local decay `(1,2)`, which corresponds to
`Dₓ⁺ → D⁰π⁺` and matches the cross-check angle in the input event.

````julia
programs = helicity_angle_programs(topology)

for (label, program) in (
    ("(((1,2),3),4)", programs[1]),
    ("((1,2),3)", programs[2]),
    ("(1,2)", programs[3]),
)
    println("program for ", label)
    foreach(step -> println("  ", format_instruction(step)), program)
end

x = cascade_kinematics(topology, system, objs)

root_angles = vertex_angles(topology, x, (((1, 2), 3), 4))
psi_angles = vertex_angles(topology, x, ((1, 2), 3))
Dx_angles = vertex_angles(topology, x, (1, 2))

println("angles for (((1,2),3),4): cosθ = ", root_angles.cosθ, ", ϕ = ", root_angles.ϕ)
println("angles for ((1,2),3): cosθ = ", psi_angles.cosθ, ", ϕ = ", psi_angles.ϕ)
println("angles for (1,2): cosθ = ", Dx_angles.cosθ, ", ϕ = ", Dx_angles.ϕ)

@assert isapprox(Dx_angles.cosθ, reference_cosθ_D_in_Dx; atol = 2e-10)
@assert isapprox(Dx_angles.ϕ, reference_ϕ_D_in_Dx; atol = 2e-10)

@assert line_invariant(topology, x, (1, 2)) ≈ mass(P_Dx)^2
@assert line_invariant(topology, x, ((1, 2), 3)) ≈ mass(P_psi)^2
@assert vertex_masses2(topology, x, (((1, 2), 3), 4)) == (mass(P_B)^2, mass(P_psi)^2, mass(pKplus)^2)
````

````
program for (((1,2),3),4)
  ToHelicityFrame((1, 2, 3, 4))
  MeasureCosThetaPhi(:v1, (1, 2, 3))
program for ((1,2),3)
  ToHelicityFrame((1, 2, 3, 4))
  ToHelicityFrame((1, 2, 3))
  MeasureCosThetaPhi(:v2, (1, 2))
program for (1,2)
  ToHelicityFrame((1, 2, 3, 4))
  ToHelicityFrame((1, 2, 3))
  ToHelicityFrame((1, 2))
  MeasureCosThetaPhi(:v3, (1,))
angles for (((1,2),3),4): cosθ = -0.04377817152173917, ϕ = 2.5569510595758183
angles for ((1,2),3): cosθ = -0.551937811945214, ϕ = 1.2182839270106544e-7
angles for (1,2): cosθ = 0.8746538492596709, ϕ = -2.849013640395373

````

## 5. Static payloads: vertices and propagators

`HadronicLineshapes.BreitWigner` is a callable propagator object. The requested
model parameters are:

```json
"Psi(4040)_mass": 4.039,
"Psi(4040)_width": 0.08
```

The left side of each pair is a bracket address, so no user-facing model input
depends on internal line or vertex numbering. The constructor resolves these
addresses once and stores typed `SVector`s internally.

For curiosity, the topology is internally compiled into a flat connectivity
matrix: rows are particle/state lines, columns are binary decay vertices, `-1`
marks the incoming parent line, and `+1` marks the two outgoing child lines.
This is the canonical graph representation from which routing and evaluation
order are derived.

````julia
connectivity = Matrix(relation(topology))
````

````
7×3 Matrix{Int64}:
  0   0   1
  0   0   1
  0   1   0
  1   0   0
  0   1  -1
  1  -1   0
 -1   0   0
````

Vertices are `ThreeBodyDecays.VertexFunction` objects. In this minimal example
we choose LS recouplings compatible with the spin assignments below.

````julia
vertices = (
    (((1, 2), 3), 4) => VertexFunction(RecouplingLS((2, 2))), # B⁺ -> ψ K
    ((1, 2), 3) => VertexFunction(RecouplingLS((0, 2))),      # ψ -> Dₓ D
    (1, 2) => VertexFunction(RecouplingLS((2, 0))),           # Dₓ -> D⁰ π
)

propagators = (
    (1, 2) => (two_j = 2, lineshape = ConstantLineshape(1.0 + 0.0im)),
    ((1, 2), 3) => (two_j = 2, lineshape = BreitWigner(4.039, 0.08)),
)

chain = DecayChain(topology; propagators, vertices);
````

## 6. General amplitude evaluation

The package-level `amplitude` method now performs the internal-helicity sum,
so the user only provides external helicities in the order `(1,2,3,4,0)`.

In formula form,

$$
\mathcal{A}(\lambda_1,\lambda_2,\lambda_3,\lambda_4,\lambda_0)
= \sum_{\lambda_{12},\lambda_{123}}
\prod_v D_v^*\,V_v \prod_r P_r.
$$

The routing contract is:

- topology chooses the local parent/children at each binary decay
- `CascadeKinematics` supplies three local masses squared
- `external_two_λs` supplies only final and initial helicities
- `CascadeSystem` and `DecayChain` supply local external/internal spins
- Julia dispatch selects vertex and propagator computations

````julia
external_two_λs = (0, 0, 0, 0, 0)
A = amplitude(chain, system, x, external_two_λs)
````

````
-0.047365002876496934 + 0.014317306220053103im
````

## 7. Numerical result

````julia
println("m²(1,2) = ", line_invariant(topology, x, (1, 2)))
println("m²((1,2),3) = ", line_invariant(topology, x, ((1, 2), 3)))
println("angles at (1,2) = ", vertex_angles(topology, x, (1, 2)))
println("summed amplitude = ", A)

@assert isfinite(real(A))
@assert isfinite(imag(A))

A
````

````
-0.047365002876496934 + 0.014317306220053103im
````

