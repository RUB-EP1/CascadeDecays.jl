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
````

````
Precompiling CascadeDecays...
   3334.8 ms  ✓ CascadeDecays
  1 dependency successfully precompiled in 4 seconds. 61 already precompiled.

````

## 1. Event four-vectors

The final-state line order is chosen to match the `InstructionalDecayTrees.jl`
angle notebook:

```text
line 1: D⁰
line 2: π⁺
line 3: D⁻
line 4: K⁺
line 5: Dₓ⁺ = D⁰π⁺
line 6: ψ(4040) = Dₓ⁺D⁻
line 7: B⁺
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
@assert vertex_lines(topology, 1) == (7, 6, 4)
@assert vertex_lines(topology, 2) == (6, 5, 3)
@assert vertex_lines(topology, 3) == (5, 1, 2)
````

## 3. Static system information

The static system stores line spins and fixed external/root masses. Internal
line masses are not static; they are supplied by the event kinematics.

We use pseudoscalar external particles, a vector `ψ(4040)`, and a vector `Dₓ`.
Spins are represented as doubled integers.

````julia
two_js = (
    0, # D⁰
    0, # π⁺
    0, # D⁻
    0, # K⁺
    2, # Dₓ⁺
    2, # ψ(4040)
    0, # B⁺
)

final_masses2 = mass.(objs) .^ 2
system = CascadeSystem(two_js, final_masses2; root_mass2 = mass(P_B)^2)
````

````
CascadeDecays.CascadeSystem{7, 4, Float64}([0, 0, 0, 0, 2, 2, 0], [3.4775909288999984, 0.019479893764752097, 3.495591122499999, 0.24371698032900024], 27.87143120480883)
````

## 4. Runtime kinematic input from the topology

`CascadeKinematics` is the completed graph-indexed input point:

- `line_masses2[line]`
- `vertex_angles[vertex]`

The angle programs are generated directly from the topology. For this chain,
the third vertex is `Dₓ⁺ → D⁰π⁺`, so it matches the cross-check angle in the
input event.

````julia
programs = helicity_angle_programs(topology)
x = cascade_kinematics(topology, system, objs)

@assert isapprox(x.vertex_angles[3].cosθ, reference_cosθ_D_in_Dx; atol = 2e-10)
@assert isapprox(x.vertex_angles[3].ϕ, reference_ϕ_D_in_Dx; atol = 2e-10)

@assert x.line_masses2[5] ≈ mass(P_Dx)^2
@assert x.line_masses2[6] ≈ mass(P_psi)^2
@assert vertex_masses2(topology, x, 1) == (x.line_masses2[7], x.line_masses2[6], x.line_masses2[4])
````

## 5. Static payloads: vertices and propagators

`HadronicLineshapes.BreitWigner` is a callable propagator object. The requested
model parameters are:

```json
"Psi(4040)_mass": 4.039,
"Psi(4040)_width": 0.08
```

Vertices are `ThreeBodyDecays.VertexFunction` objects. In this minimal example
we choose LS recouplings compatible with the static line spins above.

````julia
vertices = (
    VertexFunction(RecouplingLS((2, 2))), # B⁺ -> ψ K
    VertexFunction(RecouplingLS((0, 2))), # ψ -> Dₓ D
    VertexFunction(RecouplingLS((2, 0))), # Dₓ -> D⁰ π
)

propagators = (
    ConstantLineshape(1.0 + 0.0im),      # Dₓ⁺ line, kept non-resonant here
    BreitWigner(4.039, 0.08),            # ψ(4040) line
)

propagating_line_ids = (5, 6)

chain = DecayChain(topology, propagators, vertices, propagating_line_ids)
````

````
CascadeDecays.DecayChain{4, 2, 3, Any, ThreeBodyDecays.VertexFunction{ThreeBodyDecays.RecouplingLS, ThreeBodyDecays.NoFormFactor}, CascadeDecays.DecayTopology{7, 3, 4, Int64, 21}}(CascadeDecays.DecayTopology{7, 3, 4, Int64, 21}([0 0 1; 0 0 1; 0 1 0; 1 0 0; 0 1 -1; 1 -1 0; -1 0 0], 7, [1, 2, 3, 4]), Any[CascadeDecays.ConstantLineshape{ComplexF64}(1.0 + 0.0im), HadronicLineshapes.BreitWigner
  m: Float64 4.039
  Γ: Float64 0.08
  ma: Float64 0.0
  mb: Float64 0.0
  l: Int64 0
  d: Float64 1.0
], ThreeBodyDecays.VertexFunction{ThreeBodyDecays.RecouplingLS, ThreeBodyDecays.NoFormFactor}[ThreeBodyDecays.VertexFunction{ThreeBodyDecays.RecouplingLS, ThreeBodyDecays.NoFormFactor}(ThreeBodyDecays.RecouplingLS
  two_ls: Tuple{Int64, Int64}
, ThreeBodyDecays.NoFormFactor()), ThreeBodyDecays.VertexFunction{ThreeBodyDecays.RecouplingLS, ThreeBodyDecays.NoFormFactor}(ThreeBodyDecays.RecouplingLS
  two_ls: Tuple{Int64, Int64}
, ThreeBodyDecays.NoFormFactor()), ThreeBodyDecays.VertexFunction{ThreeBodyDecays.RecouplingLS, ThreeBodyDecays.NoFormFactor}(ThreeBodyDecays.RecouplingLS
  two_ls: Tuple{Int64, Int64}
, ThreeBodyDecays.NoFormFactor())], [5, 6])
````

## 6. General amplitude evaluation

The package-level `amplitude` method shows the intended internal contract:

- topology chooses `(parent, child1, child2)`
- `CascadeKinematics` supplies three local masses squared
- `two_λs` supplies local helicities by line id
- `CascadeSystem` supplies local spins by line id
- Julia dispatch selects vertex and propagator computations

````julia
function helicity_config(two_λψ, two_λDx)
    return (
        0,       # D⁰
        0,       # π⁺
        0,       # D⁻
        0,       # K⁺
        two_λDx,
        two_λψ,
        0,       # B⁺
    )
end

two_spin1_helicities = -2:2:2

amplitude_terms = [
    amplitude(chain, system, x, helicity_config(two_λψ, two_λDx))
    for two_λψ in two_spin1_helicities,
        two_λDx in two_spin1_helicities
]

A = sum(amplitude_terms)
````

````
-0.04736500287649317 + 0.014317306220051979im
````

## 7. Numerical result

````julia
println("line_masses2 = ", x.line_masses2)
println("vertex_angles = ", x.vertex_angles)
println("ψ(4040) BW at m²(DₓD) = ", chain.propagators[2](x.line_masses2[6]))
println("amplitude terms = ", amplitude_terms)
println("summed amplitude = ", A)

@assert isfinite(real(A))
@assert isfinite(imag(A))

A
````

````
-0.04736500287649317 + 0.014317306220051979im
````

