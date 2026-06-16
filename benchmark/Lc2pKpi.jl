# Lc2pKpi smoke script: Λc⁺ → p K⁻ π⁺ (1000 RemboOnDiet events).
# Run: julia --project=benchmark benchmark/Lc2pKpi.jl
#
# Not a performance benchmark. CascadeDecays.amplitude lacks the helicity-frame
# Wigner rotations in ThreeBodyDecays.amplitude(dc, σs); do not compare timings.

using CascadeDecays
using CascadeDecays:
    CascadeSystem,
    ConstantLineshape,
    DecayTopology,
    Propagator,
    SystemMasses,
    SystemSpinParities,
    amplitude,
    CascadeKinematics,
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

println("Lc2pKpi smoke script (shapes only — not a timing benchmark)")
println("CascadeDecays: ", pkgversion(CascadeDecays))
println("Note: Cascade amplitude omits helicity-frame Wigner rotations present in ThreeBodyDecays.amplitude(dc, σs)")
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
        propagators = (((1, 2) => Propagator(jp"1-", UNIT_LINESHAPE)),)
        minimal_ls_decay_chain(topology, system, propagators)
    end

    x = let point = events[1]
        CascadeKinematics(cascade_chain.topology, Tuple(point.momenta))
    end
    σ1 = mandelstam_lc(events[1])

    println("\n=== ", label, " ===")
    A = amplitude(cascade_chain, system, x)
    println("  Cascade external-helicity size: ", size(A), "  (p, K, π, Λc axes)")
    println("  Cascade Λc–p block A[:,1,1,:] size: ", size(A[:, 1, 1, :]))
    println("  ThreeBodyDecays amplitude size: ", size(tbd_amplitude(tbd_chain, σ1)))
end

let
    label = "Λ(1520) (pπ isobar, K spectator)"
    tbd_chain = tbd.lambda
    cascade_chain = let
        topology = DecayTopology(((1, 3), 2))
        propagators = (((1, 3) => Propagator(jp"3/2-", UNIT_LINESHAPE)),)
        minimal_ls_decay_chain(topology, system, propagators)
    end

    x = let point = events[1]
        CascadeKinematics(cascade_chain.topology, Tuple(point.momenta))
    end
    σ1 = mandelstam_lc(events[1])

    println("\n=== ", label, " ===")
    A = amplitude(cascade_chain, system, x)
    println("  Cascade external-helicity size: ", size(A), "  (p, K, π, Λc axes)")
    println("  Cascade Λc–p block A[:,1,1,:] size: ", size(A[:, 1, 1, :]))
    println("  ThreeBodyDecays amplitude size: ", size(tbd_amplitude(tbd_chain, σ1)))
end
