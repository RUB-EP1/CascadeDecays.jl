# Lc2pKpi benchmark: Λc⁺ → p K⁻ π⁺ (1000 RemboOnDiet events).
# Run: julia --project=benchmark benchmark/Lc2pKpi.jl

using BenchmarkTools
using Statistics
using CascadeDecays
using CascadeDecays:
    CascadeSystem,
    ConstantLineshape,
    DecayTopology,
    PropagatorFunction,
    SystemHelicities,
    SystemMasses,
    SystemSpinParities,
    amplitude,
    cascade_kinematics,
    minimal_ls_decay_chain
using FourVectors
using HadronicLineshapes
using Random
using RemboOnDiet
using StaticArrays
using ThreeBodyDecays:
    DecayChainLS,
    ThreeBodyMasses,
    ThreeBodyParities,
    ThreeBodySpins,
    ThreeBodySystem,
    @jp_str
import ThreeBodyDecays: amplitude as tbd_amplitude

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

# Final-state order: 1 = p, 2 = K⁻, 3 = π⁺; parent Λc⁺ (jp0 uses 1/2+).
const LC_QUANTUM = SystemSpinParities("1/2+", "0-", "0-"; jp0 = "1/2+")
const TBD_REFERENCE_TWO_ΛS = (1, 0, 0, 1)

println("Lc2pKpi amplitude benchmark")
println("CascadeDecays: ", pkgversion(CascadeDecays))
println("Dispatch: external Ne = Nf+1, internal Ni = Np (compile-time on DecayChain{Nf,Np})")
println("Events: ", N_EVENTS, " (RemboOnDiet PhaseSpaceGenerator)")
println("Quantum numbers: p 1/2+, K 0-, π 0-, Λc 1/2+ (expected external shape (2,1,1,2))")

events, tbs, Ps, tbd, system, mandelstam_lc = let
    function mandelstam_lc(point)
        p, K, π = point.momenta
        return (σ1 = mass(K + π)^2, σ2 = mass(p + π)^2, σ3 = mass(p + K)^2)
    end

    rng = MersenneTwister(0x4C633270)
    generator = PhaseSpaceGenerator([LC_MASSES.p, LC_MASSES.K, LC_MASSES.π], LC_MASSES.m0)
    events = [rand(rng, generator) for _ in 1:N_EVENTS]

    ms = ThreeBodyMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0)
    # `two_h0=...` for doubled spins; `h0=1` would x2 all entries (spin-1, not 1/2).
    spins = ThreeBodySpins(1, 0, 0; two_h0 = 1)
    tbs = ThreeBodySystem(ms, spins)
    Ps = ThreeBodyParities('+', '-', '-'; P0 = '+')

    tbd = (
        kstar = DecayChainLS(; k = 3, Xlineshape = unit_X, jp = jp"1-", Ps, tbs),
        lambda = DecayChainLS(; k = 2, Xlineshape = unit_X, jp = jp"3/2-", Ps, tbs),
    )

    masses = SystemMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0)
    system = CascadeSystem(LC_QUANTUM, masses)

    events, tbs, Ps, tbd, system, mandelstam_lc
end

let
    label = "K* (pK isobar, π spectator)"
    tbd_chain = tbd.kstar
    cascade_chain = let
        topology = DecayTopology(((1, 2), 3))
        propagators = (((1, 2) => PropagatorFunction(jp"1-", UNIT_LINESHAPE)),)
        minimal_ls_decay_chain(topology, system, propagators)
    end

    cascade = let
        xs = Vector{CascadeKinematics}(undef, length(events))
        for (i, point) in pairs(events)
            xs[i] = cascade_kinematics(cascade_chain.topology, system, Tuple(point.momenta))
        end
        hel = SystemHelicities(system.quantum.spins, 1, 0, 0; two_h0 = 1)
        (; xs, hel)
    end

    σs = mandelstam_lc.(events)

    bench_label(b) = begin
        t = median(b.times)
        string(
            round(t / 1e6; digits = 4),
            " ms total, ",
            round(t / 1e3 / N_EVENTS; digits = 3),
            " us/event",
        )
    end

    tbd_full = @benchmark begin
        for σ in $σs
            $tbd_amplitude($tbd_chain, σ)
        end
    end
    tbd_scalar = @benchmark begin
        for σ in $σs
            $tbd_amplitude($tbd_chain, σ, $TBD_REFERENCE_TWO_ΛS)
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
    A = amplitude(cascade_chain, system, cascade.xs[1])
    println("  external-helicity array size: ", size(A), "  (p, K, π, Λc axes)")
    println("  Λc–p helicity block A[:,1,1,:] size: ", size(A[:, 1, 1, :]))
    println("  ThreeBodyDecays matrix size: ", size(tbd_amplitude(tbd_chain, σs[1])))
    println("ThreeBodyDecays  full matrix : ", bench_label(tbd_full))
    println("ThreeBodyDecays  one helicity  : ", bench_label(tbd_scalar))
    println("CascadeDecays    full matrix : ", bench_label(cas_full))
    println("CascadeDecays    one helicity  : ", bench_label(cas_scalar))
    ratio = median(cas_full.times) / median(tbd_full.times)
    println("Cascade / ThreeBody (full) ≈ ", round(ratio; digits = 3), "x")
end

let
    label = "Λ(1520) (pπ isobar, K spectator)"
    tbd_chain = tbd.lambda
    cascade_chain = let
        topology = DecayTopology(((1, 3), 2))
        propagators = (((1, 3) => PropagatorFunction(jp"3/2-", UNIT_LINESHAPE)),)
        minimal_ls_decay_chain(topology, system, propagators)
    end

    cascade = let
        xs = Vector{CascadeKinematics}(undef, length(events))
        for (i, point) in pairs(events)
            xs[i] = cascade_kinematics(cascade_chain.topology, system, Tuple(point.momenta))
        end
        hel = SystemHelicities(system.quantum.spins, 1, 0, 0; two_h0 = 1)
        (; xs, hel)
    end

    σs = mandelstam_lc.(events)

    bench_label(b) = begin
        t = median(b.times)
        string(
            round(t / 1e6; digits = 4),
            " ms total, ",
            round(t / 1e3 / N_EVENTS; digits = 3),
            " us/event",
        )
    end

    tbd_full = @benchmark begin
        for σ in $σs
            $tbd_amplitude($tbd_chain, σ)
        end
    end
    tbd_scalar = @benchmark begin
        for σ in $σs
            $tbd_amplitude($tbd_chain, σ, $TBD_REFERENCE_TWO_ΛS)
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
    A = amplitude(cascade_chain, system, cascade.xs[1])
    println("  external-helicity array size: ", size(A), "  (p, K, π, Λc axes)")
    println("  Λc–p helicity block A[:,1,1,:] size: ", size(A[:, 1, 1, :]))
    println("  ThreeBodyDecays matrix size: ", size(tbd_amplitude(tbd_chain, σs[1])))
    println("ThreeBodyDecays  full matrix : ", bench_label(tbd_full))
    println("ThreeBodyDecays  one helicity  : ", bench_label(tbd_scalar))
    println("CascadeDecays    full matrix : ", bench_label(cas_full))
    println("CascadeDecays    one helicity  : ", bench_label(cas_scalar))
    ratio = median(cas_full.times) / median(tbd_full.times)
    println("Cascade / ThreeBody (full) ≈ ", round(ratio; digits = 3), "x")
end
