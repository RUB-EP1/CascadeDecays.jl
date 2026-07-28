include(joinpath(@__DIR__, "lb2lc3pi_vortices.jl"))

using FourVectors
using JSON
using LinearAlgebra
using Printf

const MR_MIN = M_LC + 2M_PI
const MR_MAX = M_LB - M_PI

const FINAL_MASSES = (M_LC, M_PI, M_PI, M_PI)
const FINAL_LABELS = ("Λc", "π⁺", "π⁻₃", "π⁻₄")

"""
    VortexChart

A topology-adapted five-dimensional chart for a sequential four-body decay,

    parent → (a,(b,c)) + d.

The three invariant-mass coordinates are `m_ab`, `m_ac`, and `m_abc`. The two
remaining coordinates are the polar angle of `a` in the `(abc)` frame and the
azimuth of `b` in the `(bc)` frame. `pair_fraction` and `cos_helicity` provide
a square chart for every fixed value of `m_abc`.
"""
struct VortexChart
    key::Symbol
    bracket::String
    bachelor::Int
    pair::NTuple{2, Int}
    spectator::Int
end

const VORTEX_CHARTS = (
    VortexChart(:a1_34, "((1,(3,4)),2)", 1, (3, 4), 2),
    VortexChart(:sigmac_123, "(((1,2),3),4)", 3, (1, 2), 4),
    VortexChart(:sigmac_124, "(((1,2),4),3)", 4, (1, 2), 3),
    VortexChart(:a1_23, "((1,(2,3)),4)", 1, (2, 3), 4),
    VortexChart(:a1_24, "((1,(2,4)),3)", 1, (2, 4), 3),
)

vortex_chart(key::Symbol) = only(filter(chart -> chart.key == key, VORTEX_CHARTS))
vortex_chart(key::AbstractString) = vortex_chart(Symbol(key))

function chart_description(chart::VortexChart)
    a, (b, c), d = chart.bachelor, chart.pair, chart.spectator
    mass_label(indices...) = "m" * join(sort(collect(indices)))
    return (
        topology = chart.bracket,
        mass_labels = (mass_label(a, b), mass_label(a, c), mass_label(a, b, c)),
        angle_labels = ("θ$(a)", "ϕ$(b)"),
        particle_labels = (
            bachelor = FINAL_LABELS[a],
            pair = (FINAL_LABELS[b], FINAL_LABELS[c]),
            spectator = FINAL_LABELS[d],
        ),
    )
end

function chart_limits(chart::VortexChart)
    a, (b, c), d = chart.bachelor, chart.pair, chart.spectator
    lower = FINAL_MASSES[a] + FINAL_MASSES[b] + FINAL_MASSES[c]
    upper = M_LB - FINAL_MASSES[d]
    return (; lower, upper)
end

"""
    chart_momenta(chart, mabc, pair_fraction, cos_helicity, theta, phi)

Construct `Λb → (a,(b,c)) d` in a topology-adapted chart. The `(abc)` system
points along +z in the Λb frame. In its rest frame particle `a` has azimuth
zero. `pair_fraction ∈ [0,1]` maps linearly onto the allowed `m_bc` interval,
and `cos_helicity` is the polar angle of particle `b` in the `(bc)` frame.

A final, fixed global rotation avoids coordinate-pole singularities in the
helicity-frame reconstruction. It changes no masses, relative angles, or zeros.
"""
function chart_momenta(
        chart::VortexChart,
        mabc,
        pair_fraction,
        cos_helicity,
        theta,
        phi,
    )
    a, (b, c), d = chart.bachelor, chart.pair, chart.spectator
    ma, mb, mc, md = FINAL_MASSES[a], FINAL_MASSES[b], FINAL_MASSES[c], FINAL_MASSES[d]
    limits = chart_limits(chart)
    limits.lower < mabc < limits.upper ||
        throw(DomainError(mabc, "cluster mass is outside phase space"))
    0 < pair_fraction < 1 ||
        throw(DomainError(pair_fraction, "pair fraction must be interior to (0,1)"))
    -1 < cos_helicity < 1 ||
        throw(DomainError(cos_helicity, "cosine must be interior to (-1,1)"))

    mbc = mb + mc + pair_fraction * (mabc - ma - mb - mc)

    q0 = breakup_momentum(M_LB, mabc, md)
    R_root = FourVector(0.0, 0.0, q0; M = mabc)
    p_d_root = FourVector(0.0, 0.0, -q0; M = md)

    qa = breakup_momentum(mabc, ma, mbc)
    p_a_R = Ry(FourVector(0.0, 0.0, qa; M = ma), theta)
    Ebc_R = sqrt(qa^2 + mbc^2)

    qb = breakup_momentum(mbc, mb, mc)
    chi = acos(cos_helicity)
    p_b_S = Rz(Ry(FourVector(0.0, 0.0, qb; M = mb), chi), phi)
    p_c_S = FourVector(-p_b_S.px, -p_b_S.py, -p_b_S.pz; M = mc)

    # Boost from the (bc) rest frame into the (abc) frame, then orient its
    # +z helicity axis along the (bc) momentum, opposite particle a.
    gamma_bc = Ebc_R / mbc
    to_R(p) = Ry(Bz(p, gamma_bc), theta + π)
    p_b_R, p_c_R = to_R(p_b_S), to_R(p_c_S)

    gamma_R = R_root.E / mabc
    to_root(p) = Bz(p, gamma_R)
    momenta = fill(to_root(p_a_R), 4)
    momenta[a] = to_root(p_a_R)
    momenta[b] = to_root(p_b_R)
    momenta[c] = to_root(p_c_R)
    momenta[d] = p_d_root

    # Numerically generic global orientation; physically equivalent to the
    # requested aligned representative.
    rotate_global(p) = Rz(Ry(p, 0.37), 0.29)
    return Tuple(map(rotate_global, momenta))
end

function chart_masses(chart::VortexChart, mabc, pair_fraction, cos_helicity)
    a, (b, c) = chart.bachelor, chart.pair
    ma, mb, mc = FINAL_MASSES[a], FINAL_MASSES[b], FINAL_MASSES[c]
    mbc = mb + mc + pair_fraction * (mabc - ma - mb - mc)
    sbc = mbc^2
    Ea = (mabc^2 - ma^2 - sbc) / (2mbc)
    Eb = (sbc + mb^2 - mc^2) / (2mbc)
    Ec = (sbc + mc^2 - mb^2) / (2mbc)
    ka = sqrt(max(0.0, Ea^2 - ma^2))
    qb = breakup_momentum(mbc, mb, mc)
    sab = ma^2 + mb^2 + 2(Ea * Eb + ka * qb * cos_helicity)
    sac = ma^2 + mc^2 + 2(Ea * Ec - ka * qb * cos_helicity)
    return (; mab = sqrt(max(0.0, sab)), mac = sqrt(max(0.0, sac)), mbc, mabc)
end

function chart_normalized_det(
        model,
        task,
        chart::VortexChart,
        mabc,
        pair_fraction,
        cos_helicity,
        theta,
        phi,
    )
    momenta = chart_momenta(chart, mabc, pair_fraction, cos_helicity, theta, phi)
    raw = amplitude(model, KinematicPoint(task, momenta))
    A = reshape(raw, 2, 2)
    return det(A) / sum(abs2, A)
end

const DEFAULT_VORTEX_CHART = vortex_chart(:a1_34)

aligned_momenta(m134, t34, coschi, theta1, phi3) =
    chart_momenta(DEFAULT_VORTEX_CHART, m134, t34, coschi, theta1, phi3)

function aligned_masses(m134, t34, coschi)
    masses = chart_masses(DEFAULT_VORTEX_CHART, m134, t34, coschi)
    return (m13 = masses.mab, m14 = masses.mac, m34 = masses.mbc)
end

aligned_normalized_det(model, task, m134, t34, coschi, theta1, phi3) =
    chart_normalized_det(
        model, task, DEFAULT_VORTEX_CHART, m134, t34, coschi, theta1, phi3,
    )

function solve_slice(f, seed; tol = 2e-9, maxiter = 45)
    y = clamp.(copy(seed), 2e-5, 1 - 2e-5)
    for _ in 1:maxiter
        fv = collect(reim(f(y)))
        norm(fv) < tol && return y, true
        J = zeros(2, 2)
        for j in 1:2
            h = min(2e-5, 0.2 * min(y[j], 1 - y[j]))
            yp, ym = copy(y), copy(y)
            yp[j] += h
            ym[j] -= h
            J[:, j] .= (collect(reim(f(yp))) .- collect(reim(f(ym)))) ./ (2h)
        end
        abs(det(J)) < 1e-9 && return y, false
        step = J \ fv
        accepted = false
        for α in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
            trial = y .- α .* step
            all(2e-5 .< trial .< 1 - 2e-5) || continue
            if abs(f(trial)) < abs(f(y))
                y = trial
                accepted = true
                break
            end
        end
        accepted || return y, false
    end
    return y, abs(f(y)) < tol
end

function roots_on_chart_slice(
        model,
        task,
        chart::VortexChart,
        mabc,
        theta,
        phi;
        seeds = Vector{Vector{Float64}}(),
        ngrid = 9,
    )
    f(y) = chart_normalized_det(
        model, task, chart, mabc, y[1], 2y[2] - 1, theta, phi,
    )
    starts = copy(seeds)
    grid = range(0.02, 0.98; length = ngrid)
    append!(starts, ([a, b] for a in grid for b in grid))
    roots = Vector{Vector{Float64}}()
    for seed in starts
        y, ok = solve_slice(f, seed)
        ok || continue
        any(r -> norm(r - y) < 2e-4, roots) || push!(roots, y)
    end
    sort!(roots; by = first)
    return roots
end

roots_on_slice(model, task, m134, theta1, phi3; kwargs...) =
    roots_on_chart_slice(
        model, task, DEFAULT_VORTEX_CHART, m134, theta1, phi3; kwargs...,
    )

function chart_dalitz_border(chart::VortexChart, mabc; n = 100)
    tvalues = range(1e-7, 1 - 1e-7; length = n)
    upper = [chart_masses(chart, mabc, t, 1 - 1e-8) for t in tvalues]
    lower = [chart_masses(chart, mabc, t, -1 + 1e-8) for t in reverse(tvalues)]
    return vcat(upper, lower, first(upper))
end

function dalitz_border(m134; n = 100)
    points = chart_dalitz_border(DEFAULT_VORTEX_CHART, m134; n)
    return [(m13 = p.mab, m14 = p.mac, m34 = p.mbc) for p in points]
end

"""
    trace_chart_strings(model, task, chart; theta, phi, nslices, ...)

Find the intersections of `det(A)=0` with a stack of fixed-`m_abc` Dalitz
planes. The full grid is used on the first plane; subsequent planes reuse the
previous roots as continuation seeds and add a smaller discovery grid. The
returned points are in the chart's three invariant masses.
"""
function trace_chart_strings(
        model,
        task,
        chart::VortexChart;
        theta = π / 3,
        phi = π / 2,
        nslices = 9,
        ngrid = 7,
        continuation_ngrid = 3,
        mlo = chart_limits(chart).lower + 2e-3,
        mhi = chart_limits(chart).upper - 2e-3,
        verbose = false,
    )
    mgrid = collect(range(mlo, mhi; length = nslices))
    slice_roots = Vector{Vector{Vector{Float64}}}(undef, nslices)
    slice_branches = Vector{Vector{Int}}(undef, nslices)
    previous = Vector{Vector{Float64}}()
    previous_branches = Int[]
    next_branch = 1
    for (i, mabc) in enumerate(mgrid)
        grid_size = i == 1 ? ngrid : continuation_ngrid
        roots = roots_on_chart_slice(
            model, task, chart, mabc, theta, phi; seeds = previous, ngrid = grid_size,
        )
        branches = fill(0, length(roots))
        candidates = sort(
            [
                (norm(previous[ip] - roots[ir]), ip, ir)
                for ip in eachindex(previous) for ir in eachindex(roots)
            ];
            by = first,
        )
        used_previous = falses(length(previous))
        for (distance, ip, ir) in candidates
            distance < 0.25 || break
            (used_previous[ip] || branches[ir] != 0) && continue
            branches[ir] = previous_branches[ip]
            used_previous[ip] = true
        end
        for ir in eachindex(branches)
            branches[ir] != 0 && continue
            branches[ir] = next_branch
            next_branch += 1
        end
        slice_roots[i] = roots
        slice_branches[i] = branches
        previous = roots
        previous_branches = branches
        verbose &&
            @printf("slice %3d/%d  mabc=%.5f  roots=%d\n", i, nslices, mabc, length(roots))
    end

    points = NamedTuple[]
    for (i, (mabc, roots, branches)) in
            enumerate(zip(mgrid, slice_roots, slice_branches))
        for (y, branch) in zip(roots, branches)
            masses = chart_masses(chart, mabc, y[1], 2y[2] - 1)
            residual = abs(chart_normalized_det(
                model, task, chart, mabc, y[1], 2y[2] - 1, theta, phi,
            ))
            push!(points, (; slice = i, branch, masses.mab, masses.mac, masses.mbc, mabc,
                pair_fraction = y[1], cos_helicity = 2y[2] - 1, residual))
        end
    end
    return (; chart, theta, phi, mgrid, slice_roots, slice_branches, points)
end

function trace_strings(
        ;
        θ1 = π / 3,
        ϕ3 = π / 2,
        nslices = 90,
        ngrid = 9,
        mlo = MR_MIN + 2e-3,
        mhi = MR_MAX - 2e-3,
    )
    model, task = lb2lc3pi_model()
    generic = trace_chart_strings(
        model,
        task,
        DEFAULT_VORTEX_CHART;
        theta = θ1,
        phi = ϕ3,
        nslices,
        ngrid,
        continuation_ngrid = ngrid,
        mlo,
        mhi,
        verbose = true,
    )
    points = [
        (; slice = p.slice, m134 = p.mabc, m13 = p.mab, m14 = p.mac, m34 = p.mbc,
            t34 = p.pair_fraction, cosχ = p.cos_helicity, p.residual)
        for p in generic.points
    ]
    return (; θ1, ϕ3, generic.mgrid, generic.slice_roots, points)
end

function main_strings()
    θ1 = deg2rad(parse(Float64, get(ENV, "VORTEX_THETA1_DEG", "60")))
    ϕ3 = deg2rad(parse(Float64, get(ENV, "VORTEX_PHI3_DEG", "90")))
    nslices = parse(Int, get(ENV, "VORTEX_NSLICES", "90"))
    ngrid = parse(Int, get(ENV, "VORTEX_NGRID", "9"))
    mlo = parse(Float64, get(ENV, "VORTEX_MR_LO", string(MR_MIN + 2e-3)))
    mhi = parse(Float64, get(ENV, "VORTEX_MR_HI", string(MR_MAX - 2e-3)))
    result = trace_strings(; θ1, ϕ3, nslices, ngrid, mlo, mhi)
    println("vortex points = ", length(result.points))
    isempty(result.points) || println("maximum residual = ", maximum(p.residual for p in result.points))
    if haskey(ENV, "VORTEX_JSON")
        payload = Dict(
            "theta1_deg" => rad2deg(result.θ1),
            "phi3_deg" => rad2deg(result.ϕ3),
            "m134_range" => [first(result.mgrid), last(result.mgrid)],
            "points" => result.points,
        )
        open(ENV["VORTEX_JSON"], "w") do io
            JSON.print(io, payload)
        end
        println("wrote ", ENV["VORTEX_JSON"])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_strings()
end
