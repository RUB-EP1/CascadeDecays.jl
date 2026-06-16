### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a1000001-0001-4000-8000-000000000002
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))

    using PlutoUI
    using FourVectors
    using ThreeBodyDecays:
        ThreeBodyMasses,
        ThreeBodyParities,
        ThreeBodySpins,
        ThreeBodySystem,
        aligned_four_vectors,
        aligned_amplitude,
        cosθij,
        Invariants,
        x2σs,
        @jp_str
    import ThreeBodyDecays: DecayChainLS, amplitude as tbd_amplitude

    include(joinpath(@__DIR__, "..", "compat", "ThreeBodyCompat.jl"))
    using .ThreeBodyCompat

    using CascadeDecays
    using CascadeDecays:
        CascadeDecay,
        CascadeSystem,
        ConstantLineshape,
        KinematicPoint,
        Propagator,
        SystemMasses,
        SystemSpinParities,
        amplitude as cascade_amplitude,
        bracket,
        kinematics_at,
        minimal_ls_decay_chain,
        vertex_angles,
        CurrentFrame,
        HelicityRootFrame
end

# ╔═╡ a1000001-0001-4000-8000-000000000001
md"""
# Λc⁺ → p K⁻ π⁺ — Λ(1520) + K\* cross-check

Compare **ThreeBodyDecays** (DecayChainLS) vs **CascadeDecays** (CascadeDecay / CascadeDecay container).

- Final states: **1 = p**, **2 = K⁻**, **3 = π⁺**
- Shared reference: **((1,2),3)** (pK isobar, π spectator)
- Λ(1520): chain **k = 3** (pK isobar)
- K\*: chain **k = 1** (Kπ isobar), evaluated with **ref_k = 3**
"""

# ╔═╡ a1000001-0001-4000-8000-000000000003
begin
const LC_MASSES = (m0 = 2.28646, p = 0.93827208816, K = 0.493677, π = 0.13957039)
const UNIT_LINESHAPE = ConstantLineshape(1.0 + 0.0im)
unit_X(σ) = UNIT_LINESHAPE(σ)

ms = ThreeBodyMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0)
Ps = ThreeBodyParities('+', '-', '-'; P0 = '+')
tbs = ThreeBodySystem(ms, ThreeBodySpins(1, 0, 0; two_h0 = 1))
system = CascadeSystem(
    SystemSpinParities("1/2+", "0-", "0-"; jp0 = "1/2+"),
    SystemMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0),
)
const REF_TOPOLOGY = three_body_topology(3)
const REF_K = three_body_k(REF_TOPOLOGY)
end

# ╔═╡ a1000001-0001-4000-8000-000000000004
md"""
## Kinematics

`aligned_four_vectors` builds the **three finals in the parent CM frame** (Σp = 0).
CascadeDecays / IDT take those three externals. For pure rotations, use `CurrentFrame()`.

We then apply the same event transform to **parent + 3 finals**:
`Rz(0.5) → Ry(0.3) → Rz(0.4)`.
"""

# ╔═╡ a1000001-0001-4000-8000-000000000014
begin
	event_transform(p) = p |> Rz(0.5) |> Ry(0.3) |> Rz(0.4)
	const EVENT_FRAME = CurrentFrame()
end

# ╔═╡ a1000001-0001-4000-8000-000000000005
@bind σ1_frac Slider(0.05:0.01:0.95, default=0.42, show_value=true)

# ╔═╡ a1000001-0001-4000-8000-000000000006
@bind σ2_frac Slider(0.05:0.01:0.95, default=0.31, show_value=true)

# ╔═╡ 536d4964-db26-4a04-9a26-3bf3eaa5b6ee
function kinematic_point(σ1_frac, σ2_frac; k = REF_K)
    σs_nt = x2σs([σ1_frac, σ2_frac], ms; k)
    finals_cm = map(
        p -> FourVector(Float64(p[1]), Float64(p[2]), Float64(p[3]); E = Float64(p[4])),
        aligned_four_vectors(σs_nt, ms; k),
    )
    p0_cm = FourVector(0.0, 0.0, 0.0; E = ms.m0)
    p0_lab = event_transform(p0_cm)
    objs = Tuple(event_transform(p) for p in finals_cm)
    σtbd = Invariants(σs_nt.σ1, σs_nt.σ2, σs_nt.σ3)
    po = reference_plane_orientation(REF_TOPOLOGY, objs; initial_frame = EVENT_FRAME)
    return (; σs_nt, σtbd, objs, po, p0_cm, p0_lab, finals_cm)
end

# ╔═╡ d32625a0-5425-473d-99ef-f3247b3fe424
global kin = kinematic_point(σ1_frac, σ2_frac)

# ╔═╡ a1000001-0001-4000-8000-000000000008
md"""
## Baseline — Λ(1520) alone
"""

# ╔═╡ a1000001-0001-4000-8000-000000000009
begin
	lambda_chain = minimal_ls_decay_chain(
	    three_body_topology(3), system,
	    ((three_body_isobar_pair(3) => Propagator(jp"3/2-", UNIT_LINESHAPE)),),
	)
	tbd_lambda = DecayChainLS(; k = 3, Xlineshape = unit_X, jp = jp"3/2-", Ps, tbs)
	baseline_task = KinematicTask(
	    (lambda_chain.topology,);
	    reference_topology = lambda_chain.topology,
	    initial_frame = EVENT_FRAME,
	)
	baseline_point = KinematicPoint(baseline_task, kin.objs)
	A_tbd_lambda_iso = tbd_amplitude(
	    tbd_lambda,
	    reference_plane_orientation(lambda_chain.topology, kin.objs; initial_frame = EVENT_FRAME),
	    kin.σtbd;
	    refζs = three_body_refζs(3),
	)
	A_cas_lambda_iso = cascade_amplitude(lambda_chain, system, baseline_point)
	(match = isapprox(A_tbd_lambda_iso, A_cas_lambda_iso), max_abs = maximum(abs, A_tbd_lambda_iso .- A_cas_lambda_iso))
end

# ╔═╡ a1000001-0001-4000-8000-00000000000a
md"""
## A. ThreeBodyDecays — DecayChainLS
"""

# ╔═╡ 145e9a3a-0673-4d3f-b122-d0f0e1fa4d64
begin
	tbd_kstar = DecayChainLS(; k = 1, Xlineshape = unit_X, jp = jp"1-", Ps, tbs)
	A_tbd_kstar = tbd_amplitude(tbd_kstar, kin.po, kin.σtbd; refζs = three_body_refζs(3))
end

# ╔═╡ a1000001-0001-4000-8000-00000000000b
A_tbd_sum = A_tbd_lambda_iso + A_tbd_kstar

# ╔═╡ a1000001-0001-4000-8000-00000000000c
md"""
## B. CascadeDecays — shared ref ((1,2),3)
"""

# ╔═╡ a1000001-0001-4000-8000-00000000000d
begin
	kstar_chain = minimal_ls_decay_chain(
	    three_body_topology(1), system,
	    ((three_body_isobar_pair(1) => Propagator(jp"1-", UNIT_LINESHAPE)),),
	)
	task = KinematicTask(
	    (lambda_chain.topology, kstar_chain.topology);
	    reference_topology = REF_TOPOLOGY,
	    wigner_finals = (1,),
	    initial_frame = EVENT_FRAME,
	)
	point = evaluate(task, kin.objs, system)
	cascade = CascadeDecay((lambda_chain, kstar_chain), system, REF_TOPOLOGY; couplings = (1.0, 1.0))
	A_cas_kstar = cascade_amplitude(kstar_chain, system, point)
end

# ╔═╡ efa2eb34-df0a-4340-871f-39b544e1e175
A_cas_sum = cascade_amplitude(cascade, point)

# ╔═╡ a1000001-0001-4000-8000-00000000000e
md"""
## Shared reference (12)3 — comparison
"""

# ╔═╡ a1000001-0001-4000-8000-00000000000f
let
function compare(label, A_tbd, A_cas)
    Δ = A_tbd .- A_cas
    (label, max_abs = maximum(abs, Δ), match = isapprox(A_tbd, A_cas))
end
    [
    compare("Λ(1520) k=3", A_tbd_lambda_iso, A_cas_lambda_iso),
    compare("K* k=1, ref_k=3", A_tbd_kstar, A_cas_kstar),
    compare("coherent sum", A_tbd_sum, A_cas_sum),
]
end

# ╔═╡ a1000001-0001-4000-8000-000000000010
md"""
## K\* diagnostics
"""

# ╔═╡ a1000001-0001-4000-8000-000000000011
kstar_diagnostics = (
    A_aligned = aligned_amplitude(tbd_kstar, kin.σtbd),
    A_tbd_zeta = tbd_amplitude(tbd_kstar, kin.σtbd; refζs = three_body_refζs(3)),
    A_tbd_po = A_tbd_kstar,
    A_cascade = A_cas_kstar,
)

# ╔═╡ a1000001-0001-4000-8000-000000000012
let d = kstar_diagnostics
    (
        max_aligned_vs_zeta = maximum(abs, d.A_aligned .- d.A_tbd_zeta),
        max_zeta_vs_cascade = maximum(abs, d.A_tbd_zeta .- d.A_cascade),
        max_po_vs_cascade = maximum(abs, d.A_tbd_po .- d.A_cascade),
    )
end

# ╔═╡ b2000001-0001-4000-8000-000000000001
md"""
## Δ diagnostics

ThreeBodyDecays `k = 2` uses the cyclic `(31)2` convention. This block keeps the
four-vectors fixed and compares the Δ chain with both `(12)3` and `(31)2`
references.
"""

# ╔═╡ b2000001-0001-4000-8000-000000000002
begin
	delta_chain = minimal_ls_decay_chain(
	    three_body_topology(2), system,
	    ((three_body_isobar_pair(2) => Propagator(jp"3/2+", UNIT_LINESHAPE)),),
	)
	tbd_delta = DecayChainLS(; k = 2, Xlineshape = unit_X, jp = jp"3/2+", Ps, tbs)

	po_ref3 = reference_plane_orientation(REF_TOPOLOGY, kin.objs; initial_frame = EVENT_FRAME)
	po_ref2 = reference_plane_orientation(delta_chain.topology, kin.objs; initial_frame = EVENT_FRAME)

	delta_point_ref3 = KinematicPoint(
	    KinematicTask(
	        (delta_chain.topology,);
	        reference_topology = REF_TOPOLOGY,
	        wigner_finals = (1,),
	        initial_frame = EVENT_FRAME,
	    ),
	    kin.objs,
	)
	delta_point_ref2 = KinematicPoint(
	    KinematicTask(
	        (delta_chain.topology,);
	        reference_topology = delta_chain.topology,
	        initial_frame = EVENT_FRAME,
	    ),
	    kin.objs,
	)

	A_tbd_delta_ref3 = tbd_amplitude(tbd_delta, po_ref3, kin.σtbd; refζs = three_body_refζs(REF_TOPOLOGY))
	A_tbd_delta_ref2 = tbd_amplitude(tbd_delta, po_ref2, kin.σtbd; refζs = three_body_refζs(delta_chain.topology))
	A_cas_delta_ref3 = cascade_amplitude(delta_chain, system, delta_point_ref3)
	A_cas_delta_ref2 = cascade_amplitude(delta_chain, system, delta_point_ref2)
end

# ╔═╡ b2000001-0001-4000-8000-000000000003
let
function compare(label, A_tbd, A_cas)
    Δ = A_tbd .- A_cas
    (label, max_abs = maximum(abs, Δ), match = isapprox(A_tbd, A_cas))
end
    [
    compare("Δ k=2, ref_k=3", A_tbd_delta_ref3, A_cas_delta_ref3),
    compare("Δ k=2, ref_k=2", A_tbd_delta_ref2, A_cas_delta_ref2),
]
end

# ╔═╡ b2000001-0001-4000-8000-000000000004
let
	x_ref2 = kinematics_at(delta_point_ref2, delta_chain.topology)
	root_angles = vertex_angles(delta_chain.topology, x_ref2, three_body_root_vertex(2))
	inner_angles = vertex_angles(delta_chain.topology, x_ref2, three_body_isobar_pair(2))
	(
	    cascade_topology = bracket(delta_chain.topology),
	    cascade_isobar = three_body_isobar_pair(2),
	    po_ref2,
	    root_vs_po = (
	        Δϕ = root_angles.ϕ - po_ref2.α,
	        Δcosθ = root_angles.cosθ - po_ref2.cosβ,
	    ),
	    inner_angles,
	)
end

# ╔═╡ b2000001-0001-4000-8000-000000000005
md"""
## Λ + K\* + Δ

For the shared `(12)3` reference, the Δ chain agrees with CascadeDecays up to an
overall sign in the ThreeBodyDecays convention. This sign comes from the
over-`2π` Wigner rotation branch, which flips spin-1/2 states. Therefore the
TBD coherent combination corresponding to the CascadeDecays `Λ + K* + Δ` model is
`Λ + K* - Δ`.
"""

# ╔═╡ b2000001-0001-4000-8000-000000000006
begin
	A_tbd_lkd_naive = A_tbd_lambda_iso + A_tbd_kstar + A_tbd_delta_ref3
	A_tbd_lkd_tbd_phase = A_tbd_lambda_iso + A_tbd_kstar - A_tbd_delta_ref3
	A_cas_lkd = cascade_amplitude(lambda_chain, system, point) + A_cas_kstar + A_cas_delta_ref3
end

# ╔═╡ b2000001-0001-4000-8000-000000000007
let
function compare(label, A_tbd, A_cas)
    Δ = A_tbd .- A_cas
    (label, max_abs = maximum(abs, Δ), match = isapprox(A_tbd, A_cas))
end
    [
    compare("naive TBD Λ + K* + Δ", A_tbd_lkd_naive, A_cas_lkd),
    compare("TBD Λ + K* - Δ", A_tbd_lkd_tbd_phase, A_cas_lkd),
]
end

# ╔═╡ a1000001-0001-4000-8000-000000000013
let
    (
    tbd_lambda = A_tbd_lambda_iso[:, 1, 1, :],
    cas_lambda = A_cas_lambda_iso[:, 1, 1, :],
    tbd_kstar = A_tbd_kstar[:, 1, 1, :],
    cas_kstar = A_cas_kstar[:, 1, 1, :],
)
end

# ╔═╡ Cell order:
# ╟─a1000001-0001-4000-8000-000000000001
# ╠═a1000001-0001-4000-8000-000000000002
# ╠═a1000001-0001-4000-8000-000000000003
# ╟─a1000001-0001-4000-8000-000000000004
# ╠═a1000001-0001-4000-8000-000000000014
# ╠═a1000001-0001-4000-8000-000000000005
# ╠═a1000001-0001-4000-8000-000000000006
# ╠═536d4964-db26-4a04-9a26-3bf3eaa5b6ee
# ╠═d32625a0-5425-473d-99ef-f3247b3fe424
# ╟─a1000001-0001-4000-8000-000000000008
# ╠═a1000001-0001-4000-8000-000000000009
# ╟─a1000001-0001-4000-8000-00000000000a
# ╠═145e9a3a-0673-4d3f-b122-d0f0e1fa4d64
# ╠═a1000001-0001-4000-8000-00000000000b
# ╟─a1000001-0001-4000-8000-00000000000c
# ╠═a1000001-0001-4000-8000-00000000000d
# ╠═efa2eb34-df0a-4340-871f-39b544e1e175
# ╟─a1000001-0001-4000-8000-00000000000e
# ╠═a1000001-0001-4000-8000-00000000000f
# ╟─a1000001-0001-4000-8000-000000000010
# ╠═a1000001-0001-4000-8000-000000000011
# ╠═a1000001-0001-4000-8000-000000000012
# ╟─b2000001-0001-4000-8000-000000000001
# ╠═b2000001-0001-4000-8000-000000000002
# ╠═b2000001-0001-4000-8000-000000000003
# ╠═b2000001-0001-4000-8000-000000000004
# ╟─b2000001-0001-4000-8000-000000000005
# ╠═b2000001-0001-4000-8000-000000000006
# ╠═b2000001-0001-4000-8000-000000000007
# ╠═a1000001-0001-4000-8000-000000000013
