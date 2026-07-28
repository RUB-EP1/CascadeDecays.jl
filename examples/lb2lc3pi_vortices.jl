using CascadeDecays
using HadronicLineshapes
using LinearAlgebra
using Random
using RamboOnDiet
using Statistics
using ThreeBodyDecays: @jp_str

const M_PI = 0.13957039
const M_LC = 2.28646
const M_LB = 5.61960

function lb2lc3pi_model()
    quantum = SystemSpinParities("1/2+", "0-", "0-", "0-"; jp0 = "1/2+")
    rho_770 = BreitWigner(0.77526, 0.1491)
    a1_1260 = BreitWigner(1.23, 0.42)
    sigmac_2455 = BreitWigner(2.45397, 20 * 0.00189)
    lcstar = BreitWigner(2.86, 5 * 0.067)

    topology_a1_23 = DecayTopology((1, ((2, 3), 4)))
    topology_a1_24 = DecayTopology((1, ((2, 4), 3)))
    topology_sigmac_123 = DecayTopology((((1, 2), 3), 4))
    topology_sigmac_124 = DecayTopology((((1, 2), 4), 3))

    chain_a1_23 = minimal_ls_decay_chain(
        topology_a1_23,
        quantum,
        (
            ((2, 3), 4) => Propagator(jp"1+", a1_1260),
            (2, 3) => Propagator(jp"1-", rho_770),
        ),
    )
    chain_a1_24 = minimal_ls_decay_chain(
        topology_a1_24,
        quantum,
        (
            ((2, 4), 3) => Propagator(jp"1+", a1_1260),
            (2, 4) => Propagator(jp"1-", rho_770),
        ),
    )
    chain_sigmac_123 = minimal_ls_decay_chain(
        topology_sigmac_123,
        quantum,
        (
            ((1, 2), 3) => Propagator(jp"3/2+", lcstar),
            (1, 2) => Propagator(jp"1/2+", sigmac_2455),
        ),
    )
    chain_sigmac_124 = minimal_ls_decay_chain(
        topology_sigmac_124,
        quantum,
        (
            ((1, 2), 4) => Propagator(jp"3/2+", lcstar),
            (1, 2) => Propagator(jp"1/2+", sigmac_2455),
        ),
    )

    model = CascadeDecay(
        (chain_a1_23, chain_a1_24, chain_sigmac_123, chain_sigmac_124),
        topology_a1_23;
        couplings = (1.0 + 0im, 1.0 + 0im, 5.0 + 0im, 5.0 + 0im),
        names = ("a1-23", "a1-24", "sigmac-123", "sigmac-124"),
    )
    task = KinematicTask(
        (topology_a1_23, topology_a1_24, topology_sigmac_123, topology_sigmac_124);
        reference_topology = topology_a1_23,
        wigner_finals = (1,),
    )
    return model, task
end

# RamboOnDiet uses eight unit-hypercube coordinates for four bodies. Removing
# the arbitrary global SO(3) orientation leaves the following five-coordinate
# section: two cluster-mass variables and three relative-angle variables.
function physical_point(u, generator)
    length(u) == 5 || throw(ArgumentError("expected five phase-space coordinates"))
    # Use a generic representative of the three redundant global rotations.
    # Putting a generated direction exactly on the z axis makes azimuthal
    # helicity frames numerically singular and creates false determinant zeros.
    r = [u[1], u[2], 0.65, 0.17, u[3], 0.43, u[4], u[5]]
    return generate_from_unit_hypercube(r, generator)
end

function amplitude_matrix(model, task, generator, u)
    event = physical_point(u, generator)
    raw = amplitude(model, KinematicPoint(task, Tuple(event.momenta)))
    size(raw) == (2, 1, 1, 1, 2) || error("unexpected amplitude shape $(size(raw))")
    return reshape(raw, 2, 2)
end

# Scale-free determinant. It has exactly the same interior zero set as det(A),
# while avoiding the large Breit-Wigner dynamic range.
function normalized_det(model, task, generator, u)
    A = amplitude_matrix(model, task, generator, u)
    return det(A) / sum(abs2, A)
end

function jacobian5(f, u; h = 2e-5)
    J = zeros(2, 5)
    for j in eachindex(u)
        hj = min(h, 0.2 * min(u[j], 1 - u[j]))
        up, um = copy(u), copy(u)
        up[j] += hj
        um[j] -= hj
        fp, fm = f(up), f(um)
        J[:, j] .= (reim(fp) .- reim(fm)) ./ (2hj)
    end
    return J
end

function solve_full(f, seed; tol = 2e-10, maxiter = 60)
    u = copy(seed)
    lo, hi = 2e-5, 1 - 2e-5
    for _ in 1:maxiter
        fv = collect(reim(f(u)))
        norm(fv) < tol && return u, true
        J = jacobian5(f, u)
        gram = J * transpose(J)
        minimum(svdvals(gram)) < 1e-14 && return u, false
        step = transpose(J) * (gram \ fv)
        accepted = false
        for α in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
            trial = u .- α .* step
            all(lo .< trial .< hi) || continue
            if abs(f(trial)) < abs(f(u))
                u = trial
                accepted = true
                break
            end
        end
        accepted || return u, false
    end
    return u, abs(f(u)) < tol
end

logit(x) = log(x / (1 - x))
logistic(x) = inv(1 + exp(-x))

function nelder_mead_minimize(g, seed; maxiter = 500, step = 0.35)
    n = length(seed)
    z0 = logit.(clamp.(seed, 1e-8, 1 - 1e-8))
    simplex = [copy(z0) for _ in 1:(n + 1)]
    for j in 1:n
        simplex[j + 1][j] += step
    end
    objective(z) = g(logistic.(z))
    vals = objective.(simplex)
    for _ in 1:maxiter
        order = sortperm(vals)
        simplex, vals = simplex[order], vals[order]
        maximum(norm(simplex[i] - simplex[1]) for i in 2:(n + 1)) < 1e-9 && break
        centroid = reduce(+, simplex[1:n]) ./ n
        reflected = centroid .+ (centroid .- simplex[end])
        fr = objective(reflected)
        if fr < vals[1]
            expanded = centroid .+ 2 .* (reflected .- centroid)
            fe = objective(expanded)
            simplex[end], vals[end] = fe < fr ? (expanded, fe) : (reflected, fr)
        elseif fr < vals[n]
            simplex[end], vals[end] = reflected, fr
        else
            contracted =
                fr < vals[end] ? centroid .+ 0.5 .* (reflected .- centroid) :
                centroid .+ 0.5 .* (simplex[end] .- centroid)
            fc = objective(contracted)
            if fc < min(fr, vals[end])
                simplex[end], vals[end] = contracted, fc
            else
                for i in 2:(n + 1)
                    simplex[i] = simplex[1] .+ 0.5 .* (simplex[i] .- simplex[1])
                    vals[i] = objective(simplex[i])
                end
            end
        end
    end
    i = argmin(vals)
    return logistic.(simplex[i]), vals[i]
end

function solve_fiber(
        f,
        fixed::NTuple{3, Float64},
        seed::NTuple{2, Float64};
        tol = 2e-10,
        maxiter = 35,
    )
    # Solve in coordinates (u4,u5), holding (u1,u2,u3) fixed.
    y = collect(seed)
    lo, hi = 2e-4, 1 - 2e-4
    for _ in 1:maxiter
        u = [fixed..., y...]
        fv = collect(reim(f(u)))
        norm(fv) < tol && return u, true
        h = 2e-5
        J = zeros(2, 2)
        for j in 1:2
            yp, ym = copy(y), copy(y)
            yp[j] += h
            ym[j] -= h
            J[:, j] .=
                (collect(reim(f([fixed..., yp...]))) .-
                 collect(reim(f([fixed..., ym...])))) ./ (2h)
        end
        abs(det(J)) < 1e-10 && return u, false
        step = J \ fv
        accepted = false
        for α in (1.0, 0.5, 0.25, 0.125, 0.0625)
            trial = y .- α .* step
            all(lo .< trial .< hi) || continue
            if abs(f([fixed..., trial...])) < abs(f(u))
                y = trial
                accepted = true
                break
            end
        end
        accepted || return u, false
    end
    u = [fixed..., y...]
    return u, abs(f(u)) < tol
end

function unique_fiber_roots(f, fixed; nseed = 9)
    roots = Vector{Vector{Float64}}()
    grid = range(0.04, 0.96; length = nseed)
    for a in grid, b in grid
        u, ok = solve_fiber(f, fixed, (a, b))
        ok || continue
        # u5 is periodic.
        distance(v, w) = hypot(v[4] - w[4], min(abs(v[5] - w[5]), 1 - abs(v[5] - w[5])))
        any(r -> distance(r, u) < 2e-4, roots) || push!(roots, u)
    end
    sort!(roots; by = u -> (u[4], u[5]))
    return roots
end

function unique_fiber_roots_nm(f, fixed; nseed = 6, tol = 2e-7)
    roots = Vector{Vector{Float64}}()
    grid = range(0.03, 0.97; length = nseed)
    fiber(y) = abs(f([fixed..., y...]))
    for a in grid, b in grid
        y, value = nelder_mead_minimize(fiber, [a, b]; maxiter = 350, step = 0.5)
        value < tol || continue
        u = [fixed..., y...]
        distance(v, w) =
            hypot(v[4] - w[4], min(abs(v[5] - w[5]), 1 - abs(v[5] - w[5])))
        any(r -> distance(r, u) < 2e-4, roots) || push!(roots, u)
    end
    sort!(roots; by = u -> (u[4], u[5]))
    return roots
end

function trace_toward(f, start, target; step = 0.015, maxiter = 600, target_tol = 2e-4)
    u = copy(start)
    path = [copy(u)]
    for _ in 1:maxiter
        distance = norm(target - u)
        distance < target_tol && return path, true
        J = jacobian5(f, u; h = 2e-6)
        tangent_projector = I - transpose(J) * ((J * transpose(J)) \ J)
        direction = tangent_projector * (target - u)
        norm(direction) < 1e-8 && return path, false
        trial = u + min(step, 0.5distance) * direction / norm(direction)
        all(2e-5 .< trial .< 1 - 2e-5) || return path, false
        corrected, ok = solve_full(f, trial; tol = 2e-9)
        ok || return path, false
        norm(target - corrected) < distance || return path, false
        u = corrected
        push!(path, copy(u))
    end
    return path, false
end

function normal_winding(f, u; radius = 1e-4, n = 128)
    V = svd(jacobian5(f, u; h = 1e-6); full = true).V
    phases = [
        angle(f(u + radius * (cos(t) * V[:, 1] + sin(t) * V[:, 2])))
            for t in range(0, 2π; length = n + 1)
    ]
    phase_change = sum(
        mod(phases[k + 1] - phases[k] + π, 2π) - π
            for k in 1:n
    )
    return round(Int, phase_change / (2π))
end

function main()
    model, task = lb2lc3pi_model()
    generator = PhaseSpaceGenerator([M_LC, M_PI, M_PI, M_PI], M_LB)
    f = u -> normalized_det(model, task, generator, u)

    rng = MersenneTwister(0x5d)
    sample = [0.001 .+ 0.998 .* rand(rng, 5) for _ in 1:5000]
    values = abs.(f.(sample))
    println("random |normalized det| quantiles = ", quantile(values, [0, 0.01, 0.1, 0.5, 0.9, 1]))

    order = sortperm(values)
    roots = Vector{Vector{Float64}}()
    for ind in order[1:12]
        u, value = nelder_mead_minimize(u -> abs(f(u)), sample[ind]; maxiter = 800)
        value < 1e-8 || continue
        any(v -> norm(u - v) < 1e-3, roots) || push!(roots, u)
    end
    println("stable interior zeros found = ", length(roots))
    for (i, u) in enumerate(roots)
        A = amplitude_matrix(model, task, generator, u)
        println(
            "  root ", i,
            ": u=", round.(u; digits = 7),
            "  |f|=", abs(f(u)),
            "  A singular-value ratio=", minimum(svdvals(A)) / maximum(svdvals(A)),
            "  J singular values=", round.(svdvals(jacobian5(f, u; h = 1e-6)); sigdigits = 6),
            "  normal winding=", normal_winding(f, u),
        )
    end

    # A nonzero pair of Jacobian singular values makes the regular zero set
    # locally 5 - 2 = 3 dimensional. Successful traces establish connectivity
    # among sampled roots; a failed target-guided trace does not prove separation.
    parent = collect(eachindex(roots))
    root_of(i) = parent[i] == i ? i : (parent[i] = root_of(parent[i]))
    function unite(i, j)
        ri, rj = root_of(i), root_of(j)
        ri == rj || (parent[rj] = ri)
    end
    for i in eachindex(roots), j in (i + 1):length(roots)
        root_of(i) == root_of(j) && continue
        connected = false
        for step in (0.01, 0.02)
            _, forward = trace_toward(f, roots[i], roots[j]; step, maxiter = 600)
            _, backward = trace_toward(f, roots[j], roots[i]; step, maxiter = 600)
            connected = forward || backward
            connected && break
        end
        connected && unite(i, j)
    end
    labels = root_of.(eachindex(roots))
    println("continuation groups among traced roots = ", length(unique(labels)), " ", labels)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
