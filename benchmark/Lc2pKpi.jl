# Lc2pKpi benchmark: Λc⁺ → p K⁻ π⁺ (1000 RemboOnDiet events).
# Run: julia --project=benchmark benchmark/Lc2pKpi.jl

using BenchmarkTools
using Statistics
using CascadeDecays
using CascadeDecays:
    CascadeSystem,
    ConstantLineshape,
    DecayChain,
    DecayTopology,
    PropagatorFunction,
    RecouplingLS,
    SystemHelicities,
    SystemMasses,
    SystemSpins,
    VertexFunction,
    amplitude,
    cascade_kinematics
using FourVectors
using HadronicLineshapes
using Random
using RemboOnDiet
using StaticArrays
using ThreeBodyDecays
using ThreeBodyDecays:
    DecayChainLS,
    MandelstamTuple,
    ThreeBodyMasses,
    ThreeBodyParities,
    ThreeBodySpins,
    ThreeBodySystem,
    amplitude as tbd_amplitude

const N_EVENTS = 1000
const UNIT_LINESHAPE = ConstantLineshape(1.0 + 0.0im)
unit_X(σ) = UNIT_LINESHAPE(σ)

# PDG-ish masses (GeV): Λc⁺, p, K⁻, π⁺
const LC_MASSES = (
    m0 = 2.28646,
    p = 0.93827208816,
    K = 0.493677,
    π = 0.13957039,
)
const KSTAR = (m = 0.89555, Γ = 0.0473, two_j = 2)   # K*(892)⁰
const L1520 = (m = 1.5195, Γ = 0.0156, two_j = 3)    # Λ(1520)⁰

function mandelstam_lc(point)
    p, K, π = point.momenta
    return (σ1 = mass(K + π)^2, σ2 = mass(p + π)^2, σ3 = mass(p + K)^2)
end

function lc_events(n::Integer, rng::AbstractRNG)
    generator = PhaseSpaceGenerator([LC_MASSES.p, LC_MASSES.K, LC_MASSES.π], LC_MASSES.m0)
    return [rand(rng, generator) for _ in 1:n]
end

function threebody_system()
    ms = ThreeBodyMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0)
    spins = ThreeBodySpins(0, 0, 0; h0 = 0)
    return ThreeBodySystem(ms, spins), ms
end

function cascade_system()
    masses = SystemMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0)
    spins = SystemSpins(0, 0, 0; two_h0 = 0)
    return CascadeSystem(spins, masses)
end

function external_helicities(system::CascadeSystem)
    return SystemHelicities(system.quantum, 0, 0, 0; two_h0 = 0)
end

"""K* in (p,K), π spectator — Cascade bracket ((1,2),3)."""
function kstar_cascade_chain()
    topology = DecayTopology(((1, 2), 3))
    return DecayChain(
        topology;
        propagators = (((1, 2) => PropagatorFunction(KSTAR.two_j, UNIT_LINESHAPE)),),
        vertices = (
            ((1, 2), 3) => VertexFunction(RecouplingLS((2, 0))),
            (1, 2) => VertexFunction(RecouplingLS((0, 2))),
        ),
    )
end

"""Λ(1520) in (p,π), K spectator — Cascade bracket ((1,3),2)."""
function lambda1520_cascade_chain()
    topology = DecayTopology(((1, 3), 2))
    return DecayChain(
        topology;
        propagators = (((1, 3) => PropagatorFunction(L1520.two_j, UNIT_LINESHAPE)),),
        vertices = (
            ((1, 3), 2) => VertexFunction(RecouplingLS((3, 0))),
            (1, 3) => VertexFunction(RecouplingLS((0, 2))),
        ),
    )
end

function threebody_reference_chains(tbs, Ps)
    kstar = DecayChainLS(;
        k = 3,
        Xlineshape = unit_X,
        jp = jp"1-",
        Ps,
        tbs,
    )
    lambda = DecayChainLS(;
        k = 2,
        Xlineshape = unit_X,
        jp = jp"3/2-",
        Ps,
        tbs,
    )
    return (; kstar, lambda)
end

function prepare_cascade(chain, system, events)
    xs = Vector{CascadeKinematics}(undef, length(events))
    for (i, point) in pairs(events)
        objs = Tuple(point.momenta)
        xs[i] = cascade_kinematics(chain.topology, system, objs)
    end
    hel = external_helicities(system)
    return (; xs, hel)
end

function bench_label(b)
    t = median(b.times)
    total_ms = t / 1e6
    per_event_us = t / 1e3 / N_EVENTS
    return string(round(total_ms; digits = 4), " ms total, ", round(per_event_us; digits = 3), " us/event")
end

function run_chain_benchmarks(; label, tbd_chain, cascade_chain, events, tbs)
    system = cascade_system()
    cascade = prepare_cascade(cascade_chain, system, events)
    σs = mandelstam_lc.(events)

    tbd_full = @benchmark begin
        for σ in $σs
            $tbd_amplitude($tbd_chain, σ)
        end
    end
    tbd_scalar = @benchmark begin
        for σ in $σs
            $tbd_amplitude($tbd_chain, σ, (0, 0, 0, 0))
        end
    end
    cas_full = @benchmark begin
        for x in $cascade.xs
            amplitude($cascade_chain, $system, x)
        end
    end
    cas_scalar = @benchmark begin
        for x in $cascade.xs
            amplitude($cascade_chain, $system, x, $cascade.hel)
        end
    end

    println("\n=== ", label, " (", N_EVENTS, " events) ===")
    println("ThreeBodyDecays  full matrix : ", bench_label(tbd_full))
    println("ThreeBodyDecays  scalar h=0  : ", bench_label(tbd_scalar))
    println("CascadeDecays    full matrix : ", bench_label(cas_full))
    println("CascadeDecays    scalar h=0  : ", bench_label(cas_scalar))
    ratio = median(cas_full.times) / median(tbd_full.times)
    println("Cascade / ThreeBody (full) ≈ ", round(ratio; digits = 3), "x")
    return (; tbd_full, tbd_scalar, cas_full, cas_scalar, ratio)
end

function main()
    println("Lc2pKpi amplitude benchmark")
    println("CascadeDecays: ", pkgversion(CascadeDecays))
    println("Dispatch: external Ne = Nf+1, internal Ni = Np (compile-time on DecayChain{Nf,Np})")
    println("Events: ", N_EVENTS, " (RemboOnDiet PhaseSpaceGenerator)")

    rng = MersenneTwister(0x4C633270)
    events = lc_events(N_EVENTS, rng)
    tbs, _ = threebody_system()
    Ps = ThreeBodyParities('+', '+', '+'; P0 = '+')
    tbd = threebody_reference_chains(tbs, Ps)

    results = Dict{String,Any}()
    results["kstar"] = run_chain_benchmarks(;
        label = "K* (pK isobar, π spectator)",
        tbd_chain = tbd.kstar,
        cascade_chain = kstar_cascade_chain(),
        events,
        tbs,
    )
    results["lambda1520"] = run_chain_benchmarks(;
        label = "Λ(1520) (pπ isobar, K spectator)",
        tbd_chain = tbd.lambda,
        cascade_chain = lambda1520_cascade_chain(),
        events,
        tbs,
    )
    return results
end

main()
