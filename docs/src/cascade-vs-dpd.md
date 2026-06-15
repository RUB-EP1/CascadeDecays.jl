```@meta
CurrentModule = CascadeDecays
```

# [CascadeDecays vs Dalitz-Plot Decomposition (TBD)](@id cascade_vs_dpd)

This note compares the `CascadeDecays.jl` helicity-frame construction with the
Dalitz-plot decomposition implemented in `ThreeBodyDecays.jl` (TBD).

The important point is that the comparison must keep the reference topology
fixed. In the examples below the reference is always

```math
((1,2),3), \qquad \mathrm{refk}=3.
```

TBD writes the three-body amplitude using the reference production-plane
orientation and Wigner small-``d`` rotations. CascadeDecays follows the full
helicity-frame paths and therefore naturally produces full Wigner
``D(\alpha,\beta,\gamma)`` rotations. For the fixed reference topology the
ordinary reference-angle phase that must be attached to the TBD amplitude is

```math
\exp\!\left(i\lambda_3\phi_{12}\right),
```

where ``\phi_{12}`` is the inner azimuth of the reference isobar ``(12)``.
This is a fixed-reference phase. It is not a target-topology spectator phase.

The identity behind the phase is

```math
D^{j*}_{\nu,\lambda_R-\lambda_3}(\phi,\theta,0)
D^{j_R*}_{\lambda_R,\lambda_1-\lambda_2}(\phi_{12},\theta_{12},0)
=
D^{j*}_{\nu,\lambda_R-\lambda_3}(\phi,\theta,\phi_{12})
\exp(i\lambda_3\phi_{12})
d^{j_R}_{\lambda_R,\lambda_1-\lambda_2}(\theta_{12}).
```

The last step comes from the rightmost Euler angle of the first ``D`` matrix,

```math
D^{j*}_{\nu,\mu}(\phi,\theta,\phi_{12})
=
D^{j*}_{\nu,\mu}(\phi,\theta,0)\,\exp(i\mu\phi_{12}),
\qquad \mu=\lambda_R-\lambda_3,
```

and from the inner isobar factor

```math
D^{j_R*}_{\lambda_R,\lambda_1-\lambda_2}(\phi_{12},\theta_{12},0)
=
\exp(-i\lambda_R\phi_{12})
d^{j_R}_{\lambda_R,\lambda_1-\lambda_2}(\theta_{12}).
```

Multiplying the two phases leaves ``\exp(i\lambda_3\phi_{12})``.

The rest of this page keeps one clean kinematic point: an aligned `(12)3`
event with physical ``\Lambda_c^+ \to p K^- \pi^+`` masses, followed by a
small but visible rotation

```julia
p |> Rz(0.1) |> Ry(0.2) |> Rz(0.3)
```

The larger angles make endpoint-``Z`` phases visible in the printed tables.

## Setup

```@example cd_vs_dpd
using LinearAlgebra
using Printf
using FourVectors
using InstructionalDecayTrees
using ThreeBodyDecays:
    ThreeBodyMasses,
    ThreeBodyParities,
    ThreeBodySpins,
    ThreeBodySystem,
    aligned_four_vectors,
    Invariants,
    x2σs,
    wr,
    cosζ,
    ispositive,
    @jp_str
import ThreeBodyDecays: DecayChainLS, amplitude as tbd_amplitude

using CascadeDecays

include(joinpath(pkgdir(CascadeDecays), "compat", "ThreeBodyCompat.jl"))
using .ThreeBodyCompat

const LC_MASSES = (m0 = 2.28646, p = 0.93827208816, K = 0.493677, π = 0.13957039)
const UNIT_LINESHAPE = ConstantLineshape(1.0 + 0.0im)
unit_X(σ) = UNIT_LINESHAPE(σ)

fourvector(p) = FourVector(Float64(p[1]), Float64(p[2]), Float64(p[3]); E = Float64(p[4]))

ms = ThreeBodyMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0)
σs_nt = x2σs([0.42, 0.31], ms; k = 3)
σs = Invariants(σs_nt.σ1, σs_nt.σ2, σs_nt.σ3)

finals_cm = fourvector.(aligned_four_vectors(σs_nt, ms; k = 3))
event_transform(p) = p |> Rz(0.1) |> Ry(0.2) |> Rz(0.3)
objs = Tuple(event_transform(p) for p in finals_cm)

const REF_K = 3
const REF_TOPOLOGY = three_body_topology(REF_K)
nothing
```

## A Reference-Phase TBD Amplitude

The next function is deliberately written in one place. It shows the three
operations needed to obtain the TBD amplitude in the convention we compare to
CascadeDecays:

1. Compute the production-plane orientation from the fixed reference topology.
2. Evaluate the TBD amplitude using `refζs` from that same reference.
3. Multiply by the fixed reference phase
   ``\exp(i\lambda_{\mathrm{ref}}\phi_{\mathrm{ref,inner}})``.

```@example cd_vs_dpd
function multiply_helicity_phase(A, axis::Integer, two_j::Integer, ϕ::Real)
    two_j == 0 && return A
    B = copy(A)
    two_λs = collect((-two_j):2:two_j)
    for I in CartesianIndices(A)
        λ = two_λs[I[axis]] / 2
        B[I] *= cis(λ * ϕ)
    end
    return B
end

function tbd_reference_amplitude(objs, tbs, chain, refk; initial_frame = CurrentFrame())
    reference_topology = three_body_topology(refk)

    # A. The reference production-plane orientation.
    po = reference_plane_orientation(reference_topology, objs; initial_frame)

    # B. TBD amplitude in the fixed reference convention.
    σ = three_body_σs(objs)
    σ_tbd = Invariants(σ.σ1, σ.σ2, σ.σ3)
    A = tbd_amplitude(chain, po, σ_tbd; refζs = three_body_refζs(refk))

    # C. The fixed-reference inner phase exp(i λ_ref ϕ_ref_inner).
    ϕ_ref_inner = po.γ
    two_j_ref = tbs.two_js[refk]
    A_phase = multiply_helicity_phase(A, refk, two_j_ref, ϕ_ref_inner)

    return (; amplitude = A_phase, raw_amplitude = A, po, ϕ_ref_inner)
end
nothing
```

For this event the reference orientation is:

```@example cd_vs_dpd
po_ref = reference_plane_orientation(REF_TOPOLOGY, objs; initial_frame = CurrentFrame())
(α = po_ref.α, cosβ = po_ref.cosβ, γ = po_ref.γ)
```

The last angle is the phase angle ``\phi_{12}``.

## External-Rotation Table

Before comparing full amplitudes, it is useful to compare the rotations alone.
The table below uses IDT path comparisons for every target topology and every
possible spin-1/2 final particle. It prints the IDT full Wigner ``D`` branch
and the TBD signed small-``d`` branch as rotation parameters only.

```@example cd_vs_dpd
idt_relative_zyz(ref_path, target_path) =
    wigner_zyz(compare_instruction_paths(ref_path, target_path, objs).relative; atol = 1e-9)
```

The TBD side has only the signed Wigner angle ``\zeta``. A positive TBD
rotation is represented as ``(0,\zeta,0)``. A negative one is represented as
``(\pi,\zeta,-\pi)``; this keeps the same small-``d`` rotation but exposes the
spinor branch.

```@example cd_vs_dpd
function tbd_relative_zyz(k, particle, σs, ms2)
    w = wr(k, REF_K, particle)
    ζ = acos(cosζ(w, σs, ms2))
    return ispositive(w) ? (ϕ = 0.0, θ = ζ, ψ = 0.0) : (ϕ = π, θ = ζ, ψ = -π)
end
```

```@example cd_vs_dpd
function path_label(path)
    join(map(path) do instr
        s = sprint(show, instr)
        kind = occursin("ToHelicityFrameParticle2", s) ? "H2" : "H1"
        m = match(r"\(\(([^)]*)\)\)", s)
        m === nothing && (m = match(r"\(([^()]*)\)", s))
        args = m === nothing ? "?" : replace(m.captures[1], ", " => "")
        "$kind($args)"
    end, " -> ")
end

function fmt_zyz(z)
    @sprintf("(% .3f, % .6f, % .3f)", z.ϕ, z.θ, z.ψ)
end

function rotation_rows()
    rows = String[]
    for k in (1, 2), particle in (1, 2, 3)
        target_topology = three_body_topology(k)
        ref_path = helicity_frame_path(REF_TOPOLOGY, particle; initial_frame = CurrentFrame())
        target_path = helicity_frame_path(target_topology, particle; initial_frame = CurrentFrame())
        idt = idt_relative_zyz(ref_path, target_path)
        tbd = tbd_relative_zyz(k, particle, σs, ms^2)

        push!(rows, @sprintf(
            "%-8s %-13d %-20s %-20s %-25s %-25s",
            k == 1 ? "(23)1" : "(31)2",
            particle,
            path_label(ref_path),
            path_label(target_path),
            fmt_zyz(idt),
            fmt_zyz(tbd),
        ))
    end
    return rows
end

println("target   spin_particle ref_path             target_path          IDT relative ZYZ          TBD relative ZYZ")
foreach(println, rotation_rows())
```

This table deliberately prints only the two rotation parameterizations. It does
not compare their spin matrices. A sign comparison belongs to the full
amplitude check below, where the fixed-reference phase
``\exp(i\lambda_3\phi_{12})`` is included.

## Full Toy Amplitude Checks

Now we repeat the exercise with simple toy amplitudes:

```math
1/2^+ \to f^+ + P^- + P^-,
```

where only one final particle is a spin-1/2 fermion at a time. The reference is
still fixed to `refk = 3`, and each target chain is compared to the
reference-phase-corrected TBD amplitude.

The isobar quantum numbers are chosen so that a chain whose spectator is the
fermion uses a vector `1-` isobar of the two pseudoscalars, while the other
chains use a fermionic `1/2-` isobar.

```@example cd_vs_dpd
function classify(A, B; atol = 1e-10)
    diff_plus = maximum(abs, A .- B)
    diff_minus = maximum(abs, A .+ B)
    sign =
        diff_plus < atol ? "+1" :
        diff_minus < atol ? "-1" :
        "not sign-only"
    return (; sign, diff_plus, diff_minus)
end

function toy_system(spin_particle)
    two_js = ntuple(i -> i == spin_particle ? 1 : 0, 3)
    tbs = ThreeBodySystem(ms, ThreeBodySpins(two_js...; two_h0 = 1))
    labels = ntuple(i -> i == spin_particle ? "1/2+" : "0-", 3)
    system = CascadeSystem(
        SystemSpinParities(labels...; jp0 = "1/2+"),
        SystemMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0),
    )
    parities = ntuple(i -> i == spin_particle ? '+' : '-', 3)
    Ps = ThreeBodyParities(parities...; P0 = '+')
    return (; two_js, tbs, system, Ps)
end

isobar_jp_for(spin_particle, k) = k == spin_particle ? jp"1-" : jp"1/2-"

function compare_toy_amplitudes(spin_particle)
    (; tbs, system, Ps) = toy_system(spin_particle)
    rows = String[]
    for k in (1, 2, 3)
        topology = three_body_topology(k)
        jp = isobar_jp_for(spin_particle, k)
        tbd_chain = DecayChainLS(; k, Xlineshape = unit_X, jp, Ps, tbs)
        cd_chain = minimal_ls_decay_chain(
            topology,
            system,
            ((three_body_isobar_pair(k) => Propagator(jp, UNIT_LINESHAPE)),),
        )
        task = KinematicTask(
            (topology,);
            reference_topology = REF_TOPOLOGY,
            wigner_finals = (1, 2, 3),
            initial_frame = CurrentFrame(),
        )
        point = KinematicPoint(task, objs)

        A_tbd = tbd_reference_amplitude(objs, tbs, tbd_chain, REF_K).amplitude
        A_cd = amplitude(cd_chain, system, point)
        c = classify(A_tbd, A_cd)
        push!(rows, @sprintf(
            "%-8d %-14s %-18.2e %-18.2e",
            k,
            c.sign,
            c.diff_plus,
            c.diff_minus,
        ))
    end
    return rows
end

for spin_particle in (1, 2, 3)
    println()
    println("spin_particle = ", spin_particle)
    println("k        sign           max |TBD-CD|       max |TBD+CD|")
    foreach(println, compare_toy_amplitudes(spin_particle))
end
```

These amplitude checks are a useful warning. The fixed reference phase
``\exp(i\lambda_3\phi_{12})`` is necessary and it fixes the particle-3 endpoint
phase in the non-reference chains. With physical LS choices the remaining
differences are global signs. The two clean central signs are:

- `(23)1` with particle 2 spin gives `-1`;
- `(31)2` with particle 1 spin gives `-1`.
