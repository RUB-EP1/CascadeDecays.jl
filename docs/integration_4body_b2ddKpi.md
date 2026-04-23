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
using InstructionalDecayTrees
using StaticArrays
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

## 2. Flat cascade topology

The incidence matrix is `relation[line, vertex]`.

- `-1`: the line enters the vertex
- `+1`: the line leaves the vertex
- `0`: the line is not attached

Vertices are ordered as:

1. `B⁺ → ψ(4040) K⁺`
2. `ψ(4040) → Dₓ⁺ D⁻`
3. `Dₓ⁺ → D⁰ π⁺`

````julia
relation = [
    0 0 1   # line 1: D⁰
    0 0 1   # line 2: π⁺
    0 1 0   # line 3: D⁻
    1 0 0   # line 4: K⁺
    0 1 -1  # line 5: Dₓ⁺
    1 -1 0  # line 6: ψ(4040)
    -1 0 0  # line 7: B⁺
]

topology = DecayTopology(relation; root = 7, finals = (1, 2, 3, 4))

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

## 4. Runtime kinematic input

`CascadeKinematics` is the completed graph-indexed input point:

- `line_masses2[line]`
- `vertex_angles[vertex]`

The final and root masses are copied from the static system. The internal
invariant masses come from this event.

````julia
program = (
    ToHelicityFrame((1, 2, 3, 4)),
    PlaneAlign((-4), (1, 2)),
    MeasureCosThetaPhi(:B_to_ψK, (1, 2, 3)),
    ToHelicityFrame((1, 2, 3)),
    MeasureCosThetaPhi(:ψ_to_DxD, (1, 2)),
    ToHelicityFrame((1, 2)),
    MeasureCosThetaPhi(:Dx_to_D0π, 1),
)

_, angle_results = apply_decay_instruction(program, objs)

@assert isapprox(angle_results.Dx_to_D0π.cosθ, reference_cosθ_D_in_Dx; atol = 2e-10)
@assert isapprox(angle_results.Dx_to_D0π.ϕ, reference_ϕ_D_in_Dx; atol = 2e-10)

internal_masses2 = (mass(P_Dx)^2, mass(P_psi)^2)
angles_by_vertex = (
    angle_results.B_to_ψK,
    angle_results.ψ_to_DxD,
    angle_results.Dx_to_D0π,
)

x = CascadeKinematics(
    topology,
    system;
    internal_masses2,
    vertex_angles = angles_by_vertex,
)

@assert x.line_masses2[5] ≈ mass(P_Dx)^2
@assert x.line_masses2[6] ≈ mass(P_psi)^2
@assert vertex_masses2(topology, x, 1) == (x.line_masses2[7], x.line_masses2[6], x.line_masses2[4])
````

## 5. Static payloads: vertices and propagators

`HadronicLineshapes.BW` is a value-level function `BW(σ, m, Γ)`, so we wrap it
in a static propagator object. The requested model parameters are:

```json
"Psi(4040)_mass": 4.039,
"Psi(4040)_width": 0.08
```

````julia
struct BreitWignerLine{T} <: AbstractLineshape
    m::T
    Γ::T
end

(p::BreitWignerLine)(σ) = HadronicLineshapes.BW(σ, p.m, p.Γ)

struct DemoVertex{T} <: AbstractVertex
    name::Symbol
    coupling::Complex{T}
end

"""
    vertex_amplitude(vertex, masses2, helicities, spins, angles)

Small dispatch-based demo vertex. It intentionally keeps the model minimal:
the vertex has a complex coupling, checks that local helicities are allowed by
the parent spin, and adds a simple helicity phase using the routed local angle.
"""
function vertex_amplitude(
    vertex::DemoVertex,
    masses2,
    helicities,
    spins,
    angles,
)
    two_λ0, two_λ1, two_λ2 = helicities
    two_j0, _, _ = spins
    two_Δλ = two_λ1 - two_λ2
    abs(two_λ0) <= two_j0 || return zero(vertex.coupling)
    abs(two_Δλ) <= two_j0 || return zero(vertex.coupling)
    return vertex.coupling * cis(0.5 * two_Δλ * angles.ϕ) * (1 + angles.cosθ) / 2
end

vertices = (
    DemoVertex(:B_to_ψK, 1.0 + 0.0im),
    DemoVertex(:ψ_to_DxD, 0.8 + 0.1im),
    DemoVertex(:Dx_to_D0π, 1.2 - 0.2im),
)

propagators = (
    ConstantLineshape(1.0 + 0.0im),      # Dₓ⁺ line, kept non-resonant here
    BreitWignerLine(4.039, 0.08),        # ψ(4040) line
)

propagating_line_ids = (5, 6)

chain = DecayChain(topology, propagators, vertices, propagating_line_ids)
````

````
CascadeDecays.DecayChain{4, 2, 3, CascadeDecays.AbstractLineshape, Main.var"##230".DemoVertex{Float64}, CascadeDecays.DecayTopology{7, 3, 4, Int64, 21}}(CascadeDecays.DecayTopology{7, 3, 4, Int64, 21}([0 0 1; 0 0 1; 0 1 0; 1 0 0; 0 1 -1; 1 -1 0; -1 0 0], 7, [1, 2, 3, 4]), CascadeDecays.AbstractLineshape[CascadeDecays.ConstantLineshape{ComplexF64}(1.0 + 0.0im), Main.var"##230".BreitWignerLine{Float64}(4.039, 0.08)], Main.var"##230".DemoVertex{Float64}[Main.var"##230".DemoVertex{Float64}(:B_to_ψK, 1.0 + 0.0im), Main.var"##230".DemoVertex{Float64}(:ψ_to_DxD, 0.8 + 0.1im), Main.var"##230".DemoVertex{Float64}(:Dx_to_D0π, 1.2 - 0.2im)], [5, 6])
````

## 6. General amplitude evaluation

This function is intentionally small and transparent. It shows the package
internal contract:

- topology chooses `(parent, child1, child2)`
- `CascadeKinematics` supplies three local masses squared
- `two_λs` supplies local helicities by line id
- `CascadeSystem` supplies local spins by line id
- Julia dispatch selects vertex and propagator computations

````julia
function routed_vertex_amplitude(chain, system, x, two_λs, vertex)
    masses2 = vertex_masses2(chain, x, vertex)
    helicities = vertex_helicities(chain, two_λs, vertex)
    spins = vertex_spins(chain.topology, system, vertex)
    angles = vertex_angles(x, vertex)
    return vertex_amplitude(chain.vertices[vertex], masses2, helicities, spins, angles)
end

function routed_propagator_product(chain, x)
    return prod(zip(chain.propagators, propagating_lines(chain))) do (propagator, line)
        propagator(line_invariant(x, line))
    end
end

function amplitude(chain, system, x, two_λs)
    V = prod(v -> routed_vertex_amplitude(chain, system, x, two_λs, v), 1:nvertices(chain))
    P = routed_propagator_product(chain, x)
    return V * P
end

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
0.174221388849182 - 0.06052072370809189im
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
0.174221388849182 - 0.06052072370809189im
````

