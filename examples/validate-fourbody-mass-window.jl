# Four-body phase-space sampling restricted to a three-pion mass window.
#
# Run from the repository root with the examples environment:
#   julia -O1 --project=examples examples/validate-fourbody-mass-window.jl [N]
# N is the event count (default 100_000, minimum 10_000). -O1 avoids an LLVM
# optimizer crash observed with Julia 1.11.5 on Apple Silicon.
#
# Setup: 190 GeV/c pi- beam on a proton at rest, outputs [recoil, pi-, pi+, pi-].
# Restricting the first RamboOnDiet cube coordinate keeps m3pi in [1.02, 1.7] GeV;
# event weights include the width of that coordinate interval. The weighted
# events estimate the restricted four-body measure delta^4(P-sum p) prod d^3p/(2E),
# without 2pi factors. They are weighted events, not an unweighted sample.
# RamboOnDiet >= v0.5.0 is required: earlier versions get the shape right but
# the absolute normalization of massive weights wrong, which this script tests.
#
# Checks: window endpoints, momentum conservation, mass shells, the raw sampling
# distribution, and the weighted m3pi spectrum against an independent two-body
# factorization integrated by Gauss-Legendre quadrature at two orders.
# Reference result (seed 20261021, N = 100_000): all assertions pass; volume
# 5.1253743 +/- 0.00422 vs quadrature 5.1311734 (1.37 sigma); effective sample size 93.65%; largest of ten mass-bin
# pulls 1.57; maximum angle discrepancy vs Cartesian boosts 6.36e-11 rad.
# Nothing is written to disk.

using RamboOnDiet, FourVectors, Random, LinearAlgebra, Printf, Test

# This is the generator's convention, without factors of 2π:
# dΦ_n = δ⁴(P-Σp) ∏ d³p/(2E). Independent, angle-integrated two-body volume.
phi2(M, a, b) = π * sqrt(max(0.0, (M-a-b)*(M+a+b)*(M-a+b)*(M+a-b))) / (2M^2)
function gauss_rule(n)
    e = eigen(SymTridiagonal(zeros(n), [k/sqrt(4k^2-1) for k in 1:n-1]))
    return e.values, 2 .* e.vectors[1,:].^2
end
function integrate(f, a, b, rule)
    x, w = rule
    return (b-a)/2 * sum(w[k]*f((a+b)/2+(b-a)/2*x[k]) for k in eachindex(x))
end
function spectrum(m, W, mp, mpi, rule)
    # dΦ4/dm = 2m Φ2(W;mp,m) Φ3(m;mpi,mpi,mpi).
    phi3 = integrate((2mpi)^2, (m-mpi)^2, rule) do t
        phi2(m, mpi, sqrt(t))*phi2(sqrt(t), mpi, mpi)
    end
    return 2m*phi2(W,mp,m)*phi3
end

# Equal pion masses; angles follow cyclic (12)3 and (23)1 conventions.
function internal_angles(m3pi, m12, m23, mpi)
    s, mu = m3pi^2, mpi^2
    sig3, sig1 = m12^2, m23^2
    sig2 = s+3mu-sig1-sig3
    lam(a,b,c) = a^2+b^2+c^2-2(a*b+a*c+b*c)
    costheta(pair, cross) = (2pair*(cross-2mu)-pair*(s-pair-mu)) /
        sqrt(max(0.0,lam(pair,mu,mu)*lam(s,pair,mu)))
    cosz = ((s+mu-sig1)*(s+mu-sig3)-2s*(sig2-2mu)) /
        sqrt(max(0.0,lam(s,mu,sig1)*lam(s,mu,sig3)))
    return acos(clamp(costheta(sig3,sig2),-1,1)),
        acos(clamp(costheta(sig1,sig3),-1,1)), -acos(clamp(cosz,-1,1))
end

# Independent Cartesian boost; these four-vectors are used only for validation.
spatial3(p) = [p.px,p.py,p.pz]
function rest_spatial(p, frame)
    beta = spatial3(frame)/frame.E
    beta2 = dot(beta,beta)
    beta2 == 0 && return spatial3(p)
    gamma = frame.E/mass(frame)
    return spatial3(p)+((gamma-1)*dot(beta,spatial3(p))/beta2-gamma*p.E)*beta
end
vector_angle(a,b) = acos(clamp(dot(a,b)/(norm(a)*norm(b)),-1,1))
function cartesian_angles(p1,p2,p3)
    total = p1+p2+p3
    return vector_angle(rest_spatial(p1,p1+p2),-rest_spatial(p3,p1+p2)),
        vector_angle(rest_spatial(p2,p2+p3),-rest_spatial(p1,p2+p3)),
        -vector_angle(rest_spatial(p1,total),rest_spatial(p3,total))
end

function validate_window(; n=100_000, seed=20261021)
    n >= 10_000 || error("Use at least 10_000 events for the histogram checks")
    mpi, mp, beam_p = 0.13957039, 0.93827208816, 190.0
    W = sqrt(mp^2 + mpi^2 + 2mp*sqrt(beam_p^2+mpi^2))
    ma, mb = 1.02, 1.7
    delta = W-mp-3mpi
    coordinate(m) = begin
        z = ((m-3mpi)/delta)^2
        z^2*(3-2z)
    end
    ua, ub = coordinate(ma), coordinate(mb)
    width = ub-ua
    gen = PhaseSpaceGenerator([mp,mpi,mpi,mpi],W)
    rng = MersenneTwister(seed)
    edges = collect(range(ma,mb;length=11))
    counts, sums, sums2 = zeros(Int,10), zeros(10), zeros(10)
    minmass, maxmass = Inf, -Inf
    shell_error = conservation_error = map_error = 0.0
    positive_weights = finite_angles = true
    angle_error = 0.0
    for _ in 1:n
        u = rand(rng,8)
        u[1] = ua+width*u[1]
        # Mass coordinates must be in the open unit interval.
        while iszero(u[2]); u[2]=rand(rng); end
        event = generate_from_unit_hypercube(u,gen)
        # Output ordering is recoil, pi-, pi+, pi-; retain this for later amplitudes.
        ps = event.momenta
        parent = RamboOnDiet.total_momentum(ps)
        m = mass(ps[2]+ps[3]+ps[4])
        weight = width*phase_space_weight(event)
        m12, m23 = mass(ps[2]+ps[3]), mass(ps[3]+ps[4])
        angles = internal_angles(m,m12,m23,mpi)
        independent = cartesian_angles(ps[2],ps[3],ps[4])
        angle_error = max(angle_error,maximum(abs.(angles .- independent)))
        finite_angles &= all(isfinite,angles)
        positive_weights &= isfinite(weight) && weight>0
        minmass, maxmass = min(minmass,m), max(maxmass,m)
        map_error = max(map_error,abs(coordinate(m)-u[1]))
        conservation_error = max(conservation_error,abs(parent.E-W),abs(parent.px),abs(parent.py),abs(parent.pz))
        shell_error = max(shell_error,maximum(abs(mass(p)-mi) for (p,mi) in zip(ps,[mp,mpi,mpi,mpi])))
        bin = clamp(searchsortedlast(edges,m),1,10)
        counts[bin] += 1
        sums[bin] += weight
        sums2[bin] += weight^2
    end
    # Compare absolute integrals, not just a histogram normalized to itself.
    reference(rule) = [integrate(m->spectrum(m,W,mp,mpi,rule),edges[b],edges[b+1],rule) for b in 1:10]
    coarse, exact = reference(gauss_rule(64)), reference(gauss_rule(128))
    estimates = sums/n
    errors = sqrt.((sums2/n-estimates.^2)/(n-1))
    pulls = (estimates-exact)./errors
    total_mc, total_ref = sum(estimates), sum(exact)
    total_se = sqrt((sum(sums2)/n-total_mc^2)/(n-1))
    # Raw (unweighted) events follow the restricted cube map, not exact massive Φ4.
    raw_prob = diff(coordinate.(edges))/width
    raw_pulls = (counts-n*raw_prob)./sqrt.(n*raw_prob.*(1 .- raw_prob))
    neff = sum(sums)^2/sum(sums2)
    @testset "Restricted four-body mass window" begin
        @test ma <= minmass <= maxmass <= mb
        @test positive_weights
        @test finite_angles
        @test angle_error < 2e-8
        @test conservation_error < 1e-9
        @test shell_error < 1e-9
        @test map_error < 1e-11
        @test maximum(abs.(coarse./exact .- 1)) < 1e-5
        @test maximum(abs.(raw_pulls)) < 5.5
        @test maximum(abs.(pulls)) < 5.5
        @test abs(total_mc-total_ref) < 5.5total_se
        # Check endpoints with different downstream coordinates.
        for m in (ma, (ma+mb)/2, mb), r in (0.1,0.5,0.9)
            u=fill(r,8); u[1]=coordinate(m)
            ps=generate_from_unit_hypercube(u,gen).momenta
            @test mass(ps[2]+ps[3]+ps[4]) ≈ m atol=1e-9
        end
    end
    @printf("Max invariant vs Cartesian angle error: %.3g rad\n",angle_error)
    println("RamboOnDiet source: ",pkgdir(RamboOnDiet))
    @printf("seed=%d, N=%d, W=%.12f GeV\n",seed,n,W)
    @printf("u1 window=[%.12g, %.12g], width=%.12g\n",ua,ub,width)
    @printf("m3pi range=[%.9f, %.9f] GeV; all %d points in window\n",minmass,maxmass,n)
    @printf("Unrestricted acceptance=%.6f%%; avoided rejection factor=%.1f\n",100width,1/width)
    @printf("Max errors: conservation=%.3g GeV, mass shell=%.3g GeV, cube map=%.3g\n",conservation_error,shell_error,map_error)
    @printf("Restricted Φ4: MC=%.8g +/- %.3g; quadrature=%.8g; pull=%.3f\n",total_mc,total_se,total_ref,(total_mc-total_ref)/total_se)
    @printf("Effective weighted sample=%.1f (%.2f%%); quadrature refinement=%.3g\n",neff,100neff/n,maximum(abs.(coarse./exact .- 1)))
    println("m3pi bin     count    weighted fraction  reference fraction   pull")
    for b in 1:10
        @printf("%.3f--%.3f  %7d       %.6f           %.6f       % .3f\n",edges[b],edges[b+1],counts[b],estimates[b]/total_mc,exact[b]/total_ref,pulls[b])
    end
    return (;total_mc,total_ref,total_se,neff,pulls,angle_error)
end

validate_window(; n=isempty(ARGS) ? 100_000 : parse(Int,ARGS[1]));
