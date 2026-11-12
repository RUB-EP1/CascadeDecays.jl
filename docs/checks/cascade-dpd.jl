# Executed by cascade-vs-dpd.qmd; also runnable with julia --project=docs.
using Random, Printf
using FourVectors, CascadeDecays
import ThreeBodyDecays as TBD
include(joinpath(@__DIR__, "..", "..", "compat", "ThreeBodyCompat.jl"))
using .ThreeBodyCompat

oriented_tbd(chain, po, invariants, r) = TBD.amplitude(chain, po, invariants; reference_k = r)
tbd_vertex(ls) = TBD.Vertex(TBD.RecouplingLS(ls))

conversion_sign(two_js, k, r) = k == r ? 1 : (-1)^two_js[6 - k - r]

function reference_phase(A, two_js, r, γ)
    B = copy(A)
    for I in CartesianIndices(B)
        λr = (-two_js[r] + 2(I[r] - 1)) / 2
        B[I] *= cis(λr * γ)
    end
    return B
end

# Enumerate CD's allowed LS pairs and give those exact pairs to TBD.
# Do not rely on the two libraries sorting their LS lists identically.
function matched_chains(two_js, two_j0, two_jR, k, ms)
    jp(n) = string(n, "/2+")
    quantum = SystemSpinParities(jp.(two_js)...; jp0 = jp(two_j0))
    topology = three_body_topology(k)
    pair = three_body_isobar_pair(k)
    propagators = (pair => Propagator(CascadeDecays.SpinParity(jp(two_jR)), ConstantLineshape(1.0 + 0im)),)
    specs = CascadeDecays.possible_vertex_couplings(topology, quantum, propagators)
    tbs = TBD.ThreeBodySystem(ms, TBD.ThreeBodySpins(two_js...; two_h0 = two_j0))
    return map(collect(Iterators.product(last(specs[1]), last(specs[2])))) do (LS, ls)
        cd = CascadeDecays.DecayChain(
            topology, quantum.spins;
            propagators,
            vertices = (
                first(specs[1]) => CascadeDecays.Vertex(CascadeDecays.RecouplingLS(LS)),
                first(specs[2]) => CascadeDecays.Vertex(CascadeDecays.RecouplingLS(ls)),
            )
        )
        tbd = TBD.DecayChain(;
            k, two_j = two_jR, tbs, Xlineshape = _ -> 1.0 + 0im,
            HRk = tbd_vertex(LS), Hij = tbd_vertex(ls)
        )
        (; cd, tbd, LS, ls)
    end
end

function sampled_event(rng, ms)
    # Interior Dalitz points, with independent, isotropic overall orientations.
    σ = TBD.x2σs(0.05 .+ 0.9 .* rand(rng, 2), ms; k = 3)
    α, γ = 2π .* rand(rng, 2) .- π
    β = acos(2rand(rng) - 1)
    return Tuple(
        begin
            p = FourVector(Float64(v[1]), Float64(v[2]), Float64(v[3]); E = Float64(v[4]))
            p |> Rz(γ) |> Ry(β) |> Rz(α)
        end for v in TBD.aligned_four_vectors(σ, ms; k = 3)
    )
end

"""
    validate_conversion(; nsamples=8, seed=20261112, cases=..., atol=2e-10, rtol=2e-10)

Compare every complex helicity component for all nine (chain, reference) pairs,
all allowed LS pairs for two isobar-spin choices, and reproducible four-vectors.
Cases contain `(final doubled spins, parent doubled spin)` and must obey spin parity.
This is a finite numerical stress test, not a proof for unbounded spins.
"""
function validate_conversion(;
        nsamples = 8, seed = 20261112,
        cases = [
            ((0, 0, 0), 0), ((1, 0, 0), 1), ((0, 1, 0), 1), ((0, 0, 1), 1),
            ((1, 1, 0), 0), ((1, 0, 1), 2), ((0, 1, 1), 2), ((1, 1, 1), 1),
            ((2, 0, 0), 2), ((3, 1, 0), 2), ((1, 2, 1), 2),
        ],
        atol = 2.0e-10, rtol = 2.0e-10
    )
    nsamples > 0 || throw(ArgumentError("nsamples must be positive"))
    rng = MersenneTwister(seed)
    ms = TBD.ThreeBodyMasses(0.93827208816, 0.493677, 0.13957039; m0 = 2.28646)
    events = [sampled_event(rng, ms) for _ in 1:nsamples]
    checks = 0; max_error = 0.0; phase_witness = false; sign_witness = false
    covered = Set{Tuple{Int, Int}}()
    for (two_js, two_j0) in cases
        iseven(two_j0 - sum(two_js)) || throw(ArgumentError("incompatible spin parity"))
        for k in 1:3
            i, j = three_body_isobar_pair(k)
            for two_jR in (two_js[i] + two_js[j], two_js[i] + two_js[j] + 2)
                chains = matched_chains(two_js, two_j0, two_jR, k, ms)
                isempty(chains) && continue # forbidden parity/triangle assignment
                for r in 1:3
                    task = KinematicTask(
                        (three_body_topology(k),);
                        reference_topology = three_body_topology(r),
                        wigner_finals = (1, 2, 3), initial_frame = CurrentFrame()
                    )
                    for (event_index, objs) in enumerate(events)
                        point = KinematicPoint(task, objs)
                        σ = three_body_σs(objs)
                        invariants = TBD.Invariants(σ.σ1, σ.σ2, σ.σ3)
                        po = reference_plane_orientation(
                            three_body_topology(r), objs;
                            initial_frame = CurrentFrame()
                        )
                        for chain in chains
                            raw = oriented_tbd(chain.tbd, po, invariants, r)
                            phased = reference_phase(raw, two_js, r, po.γ)
                            sign = conversion_sign(two_js, k, r)
                            predicted = sign .* phased
                            actual = CascadeDecays.amplitude(chain.cd, point)
                            err = maximum(abs, actual - predicted)
                            scale = max(maximum(abs, actual), maximum(abs, predicted))
                            err <= atol + rtol * scale || error(
                                "CD–TBD mismatch: seed=$seed event=$event_index spins=$two_js " *
                                    "parent=$two_j0 isobar=$two_jR k=$k r=$r " *
                                    "LS=$(chain.LS) ls=$(chain.ls) error=$err scale=$scale"
                            )
                            max_error = max(max_error, err)
                            checks += 1; push!(covered, (k, r))
                            phase_witness |= maximum(abs, actual - sign .* raw) > 1.0e-5
                            sign_witness |= sign == -1 && maximum(abs, actual - phased) > 1.0e-5
                        end
                    end
                end
            end
        end
    end
    length(covered) == 9 || error("not all chain/reference pairs were exercised")
    phase_witness || error("sample did not detect omission of the continuous phase")
    sign_witness || error("sample did not detect omission of the discrete sign")
    return (;
        seed, samples = nsamples, spin_cases = length(cases), checks, max_error,
        phase_witness, sign_witness,
    )
end

function print_validation(result)
    @printf(
        "%d tensor comparisons; %d spin assignments; %d sampled events; all 9 (k,r) pairs\n",
        result.checks, result.spin_cases, result.samples
    )
    @printf("Maximum absolute component error: %.3e\n", result.max_error)
    return println("Omitting either the phase or sign is detected. Seed: ", result.seed)
end

if abspath(PROGRAM_FILE) == @__FILE__
    print_validation(validate_conversion())
end
