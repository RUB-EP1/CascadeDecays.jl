# Run from the repository root: julia --project=. docs/validate-three-pion-dpd.jl
using CascadeDecays, FourVectors, InstructionalDecayTrees, PartialWaveFunctions
using LinearAlgebra, Random, Printf, Test
import ThreeBodyDecays as TBD

const PION_MASS = 0.13957039
const PARENT_MASS = 1.4
const MASSES_3PI = TBD.ThreeBodyMasses(PION_MASS, PION_MASS, PION_MASS; m0 = PARENT_MASS)

# Random physical events at fixed parent mass, not a uniform phase-space sampler.
function random_three_pions(rng)
    v1, v2 = randn(rng, 3), randn(rng, 3)
    vs = (v1, v2, -v1 - v2)
    energy(a) = sum(sqrt(PION_MASS^2 + a^2 * dot(v, v)) for v in vs)
    lo, hi = 0.0, 1.0
    while energy(hi) < PARENT_MASS
        hi *= 2
    end
    for _ in 1:100
        a = (lo + hi) / 2
        energy(a) > PARENT_MASS ? (hi = a) : (lo = a)
    end
    return Tuple(FourVector((lo * v)...; E = sqrt(PION_MASS^2 + lo^2 * dot(v, v))) for v in vs)
end

function chain_angles(tree, ps)
    programs = CascadeDecays.helicity_angle_programs(DecayTopology(tree); initial_frame = CurrentFrame())
    return map(programs) do program
        _, result = apply_decay_instruction(program, ps)
        a = only(values(result))
        (phi = a.ϕ, theta = acos(clamp(a.cosθ, -1, 1)))
    end
end

rz3(x) = [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1]
ry3(x) = [cos(x) 0 sin(x); 0 1 0; -sin(x) 0 cos(x)]
orientation(a) = rz3(a[1].phi) * ry3(a[1].theta) * rz3(a[2].phi)
spatial(p) = [p.px, p.py, p.pz]
parity(n) = iseven(n) ? 1.0 : -1.0
d(j, m, n, theta) = wignerd(j, m, n, cos(theta))
Dc(j, m, n, phi, theta, psi = 0.0) = conj(wignerD(j, m, n, phi, cos(theta), psi))

# Independent Cartesian boost, followed by projection onto the helicity axes.
function cartesian_angles(ps, i, j)
    pair = ps[i] + ps[j]
    q, v = spatial(pair), spatial(ps[i])
    z = normalize(q)
    phi, theta = atan(z[2], z[1]), acos(z[3])
    axes = rz3(phi) * ry3(theta)
    beta = q / pair.E
    gamma = pair.E / mass(pair)
    vstar = v + ((gamma - 1) * dot(beta, v) / dot(beta, beta) - gamma * ps[i].E) * beta
    local_v = transpose(axes) * vstar
    return (
        (phi = phi, theta = theta),
        (phi = atan(local_v[2], local_v[1]), theta = acos(local_v[3] / norm(local_v))),
    )
end

function validate_three_pion_dpd(; nevents = 100)
    rng = MersenneTwister(20261021)
    errors = Dict(
        k => 0.0 for k in (
                :angles, :zeta, :rotation13, :rotation31,
                :pull13, :pull31, :reversal, :reflectivity, :reflectivity_M0, :missing_phase, :mixed_TY,
            )
    )
    record!(key, value) = (errors[key] = max(errors[key], abs(value)))
    first_event = nothing
    @testset "Three-pion angular factorization" begin
        for event in 1:nevents
            ps = random_three_pions(rng)
            a = chain_angles(((1, 2), 3), ps)
            b = chain_angles(((1, 3), 2), ps)
            c = chain_angles(((3, 1), 2), ps)
            @test norm(spatial(sum(ps))) < 1.0e-14
            @test sum(p.E for p in ps) ≈ PARENT_MASS
            @test all(p -> isapprox(mass(p), PION_MASS; atol = 1.0e-14), ps)
            for (angles, ij) in ((a, (1, 2)), (b, (1, 3)), (c, (3, 1)))
                independent = cartesian_angles(ps, ij...)
                for v in 1:2
                    record!(:angles, angles[v].theta - independent[v].theta)
                    record!(:angles, angle(cis(angles[v].phi - independent[v].phi)))
                end
            end
            sig = (mass(ps[2] + ps[3])^2, mass(ps[1] + ps[3])^2, mass(ps[1] + ps[2])^2)
            w = TBD.wr(2, 3, 0)
            zeta = (TBD.ispositive(w) ? 1 : -1) * acos(clamp(TBD.cosζ(w, sig, MASSES_3PI^2), -1, 1))
            q = transpose(orientation(a)) * orientation(c)
            record!(:zeta, zeta - atan(q[1, 3], q[3, 3]))
            record!(:rotation13, norm(transpose(orientation(a)) * orientation(b) - ry3(zeta) * rz3(pi)))
            record!(:rotation31, norm(q - ry3(zeta)))
            @test cos(c[2].theta) ≈ TBD.cosθij(sig, MASSES_3PI^2; k = 2) atol = 1.0e-12
            @test b[2].theta + c[2].theta ≈ pi
            if event == 1
                first_event = (; ps, a, b, c, zeta, q)
            end
            for J in 0:4, S in 0:3, M in -J:J, lam in -min(J, S):min(J, S)
                direct13 = Dc(J, M, lam, b[1].phi, b[1].theta) * Dc(S, lam, 0, b[2].phi, b[2].theta)
                direct31 = Dc(J, M, lam, c[1].phi, c[1].theta) * Dc(S, lam, 0, c[2].phi, c[2].theta)
                pulled = sum(Dc(J, M, lp, a[1].phi, a[1].theta, a[2].phi) * d(J, lp, lam, zeta) for lp in -J:J)
                record!(:pull13, direct13 - pulled * parity(lam) * d(S, lam, 0, b[2].theta))
                record!(:pull31, direct13 - parity(S) * pulled * d(S, lam, 0, c[2].theta))
                record!(:reversal, direct13 - parity(S) * direct31)
                # Negative controls: omit the phase, or mix chain-12 TY with chain-13 GJ.
                record!(:missing_phase, direct13 - pulled * d(S, lam, 0, b[2].theta))
                mixed = Dc(J, M, lam, b[1].phi, b[1].theta) * Dc(S, lam, 0, a[2].phi, a[2].theta)
                record!(:mixed_TY, mixed - pulled * parity(lam) * d(S, lam, 0, b[2].theta))
            end
            # Test the full coherent sum with complex couplings. Arbitrary H values
            # test the linear basis change, independently of dynamical parity constraints.
            for J in 0:4, S in 0:3
                helicities = -min(J, S):min(J, S)
                H12(lam) = complex(0.7 + 0.1 * lam, 0.2 - 0.03 * lam)
                H13(lam) = complex(-0.4 + 0.07 * lam, 0.5 + 0.04 * lam)
                h12, h13 = 0.8 - 0.3im, -0.2 + 0.6im
                original(M) = sum(
                    Dc(J, M, lam, a[1].phi, a[1].theta) * H12(lam) * Dc(S, lam, 0, a[2].phi, a[2].theta) * h12 +
                        Dc(J, M, lam, b[1].phi, b[1].theta) * H13(lam) * Dc(S, lam, 0, b[2].phi, b[2].theta) * h13
                        for lam in helicities
                )
                O = [
                    (abs(lp) <= S ? H12(lp) * d(S, lp, 0, a[2].theta) * h12 : 0im) +
                        parity(S) * sum(d(J, lp, lam, zeta) * H13(lam) * d(S, lam, 0, c[2].theta) * h13 for lam in helicities)
                        for lp in -J:J
                ]
                for P in (-1, 1), eps in (-1, 1), M in 0:J
                    normM = M == 0 ? 0.5 : inv(sqrt(2.0))
                    coefficient = eps * P * parity(J - M)
                    direct = normM * (original(M) - coefficient * original(-M))
                    pulled = sum(
                        normM * (
                            Dc(J, M, lp, a[1].phi, a[1].theta, a[2].phi) -
                                coefficient * Dc(J, -M, lp, a[1].phi, a[1].theta, a[2].phi)
                        ) * O[lp + J + 1] for lp in -J:J
                    )
                    record!(:reflectivity, direct - pulled)
                    if M == 0
                        expected = eps == -P * parity(J) ? original(0) : 0im
                        record!(:reflectivity_M0, pulled - expected)
                    end
                end
            end
        end
        for key in (:angles, :zeta, :rotation13, :rotation31, :pull13, :pull31, :reversal, :reflectivity, :reflectivity_M0)
            @test errors[key] < 1.0e-10
        end
        @test errors[:missing_phase] > 0.1
        @test errors[:mixed_TY] > 0.1
    end
    println("Seed = 20261021; events = $nevents; J = 0:4; S = 0:3")
    println("First event four-vectors (E, px, py, pz), GeV:")
    for (i, p) in enumerate(first_event.ps)
        @printf("p%d = (% .12f, % .12f, % .12f, % .12f)\n", i, p.E, p.px, p.py, p.pz)
    end
    println("chain     phi_GJ         theta_GJ       phi_TY         theta_TY (radians)")
    for (name, a) in (("(12)3", first_event.a), ("(13)2", first_event.b), ("(31)2", first_event.c))
        @printf("%s  % .12f % .12f % .12f % .12f\n", name, a[1].phi, a[1].theta, a[2].phi, a[2].theta)
    end
    @printf("zeta = %.12f rad = %.12f deg\n", first_event.zeta, rad2deg(first_event.zeta))
    for key in sort(collect(keys(errors)); by = string)
        @printf("max %-15s %.6e\n", key, errors[key])
    end
    return (; first_event, errors)
end

validate_three_pion_dpd();
