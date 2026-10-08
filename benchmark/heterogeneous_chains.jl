# Evaluation-speed benchmark for a cascade whose chains mix payload types.
# Run: julia --project=benchmark benchmark/heterogeneous_chains.jl
#
# pp -> p p K+ K- with eight chains on four topologies. Lineshapes mix
# BreitWigner (l = 0..3), Flatte, and ConstantLineshape; vertices carry
# BlattWeisskopf{L} form factors with L from the minimal LS coupling. Such a
# chain cannot be stored in a concretely typed `SVector`, which is what this
# benchmark is meant to guard.
#
# The timing includes the HadronicLineshapes lineshape calls; `BreitWigner`
# there is itself type-unstable (~1 μs per call), so it is a sizeable part of
# the total.

using BenchmarkTools
using CascadeDecays
using HadronicLineshapes
using Random
using RamboOnDiet
using ThreeBodyDecays: @jp_str, RecouplingLS

const mp, mK = 0.9382720813, 0.493677
const quantum = SystemSpinParities("1/2+", "1/2+", "0-", "0-"; jp0 = "0+")

const_ls = ConstantLineshape(1.0 + 0.0im)
bw(m, Γ, ma, mb, l) = BreitWigner(; m, Γ, ma, mb, l, d = 1.5)
flatte_a0 = Flatte(; m = 0.98, gsq1 = 0.2, ma1 = mK, mb1 = mK, gsq2 = 0.1, ma2 = 0.14, mb2 = 0.55)

function chain_with_form_factors(topology, propagators)
    specs = minimal_vertex_couplings(topology, quantum, propagators)
    vertices = map(specs) do (address, (two_l, two_s))
        address => Vertex(RecouplingLS((two_l, two_s)), BlattWeisskopf{div(two_l, 2)}(1.5))
    end
    return DecayChain(topology, quantum.spins; propagators, vertices)
end

t_KK = DecayTopology(((1, 2), (3, 4)))
t_L1 = DecayTopology(((1, 4), (2, 3)))
t_L2 = DecayTopology(((2, 4), (1, 3)))
t_seq = DecayTopology((((1, 4), 3), 2))

chains = (
    chain_with_form_factors(t_KK, ((1, 2) => Propagator(jp"0+", const_ls), (3, 4) => Propagator(jp"1-", bw(1.019, 0.0042, mK, mK, 1)))),
    chain_with_form_factors(t_KK, ((1, 2) => Propagator(jp"0+", const_ls), (3, 4) => Propagator(jp"2+", bw(1.525, 0.073, mK, mK, 2)))),
    chain_with_form_factors(t_KK, ((1, 2) => Propagator(jp"1-", const_ls), (3, 4) => Propagator(jp"0+", flatte_a0))),
    chain_with_form_factors(t_L1, ((1, 4) => Propagator(jp"3/2-", bw(1.5195, 0.0156, mp, mK, 2)), (2, 3) => Propagator(jp"1/2-", const_ls))),
    chain_with_form_factors(t_L2, ((2, 4) => Propagator(jp"3/2-", bw(1.5195, 0.0156, mp, mK, 2)), (1, 3) => Propagator(jp"1/2-", const_ls))),
    chain_with_form_factors(t_L1, ((1, 4) => Propagator(jp"1/2-", bw(1.6, 0.15, mp, mK, 0)), (2, 3) => Propagator(jp"1/2+", const_ls))),
    chain_with_form_factors(t_L1, ((1, 4) => Propagator(jp"5/2+", bw(1.82, 0.08, mp, mK, 3)), (2, 3) => Propagator(jp"1/2-", const_ls))),
    chain_with_form_factors(
        t_seq, (
            (1, 4) => Propagator(jp"3/2-", bw(1.5195, 0.0156, mp, mK, 2)),
            ((1, 4), 3) => Propagator(jp"1/2+", bw(2.1, 0.2, 1.5195, mK, 1)),
        )
    ),
)

model = CascadeDecay(chains, t_KK; couplings = ntuple(i -> complex(1.0 / i, 0.3), length(chains)))

task = KinematicTask((t_KK, t_L1, t_L2, t_seq); reference_topology = t_KK, wigner_finals = (1, 2))
rng = MersenneTwister(1)
generator = PhaseSpaceGenerator([mp, mp, mK, mK], 3.5)
points = [KinematicPoint(task, Tuple(rand(rng, generator).momenta)) for _ in 1:200]
point = first(points)

println("CascadeDecays ", pkgversion(CascadeDecays))
println("payload storage: ", nameof(typeof(first(chains).propagators)))
for (i, chain) in enumerate(chains)
    trial = @benchmark amplitude($chain, $point)
    println("chain $i: ", round(median(trial).time / 1.0e3; digits = 2), " μs, ", trial.memory, " bytes")
end
trial = @benchmark unpolarized_intensity($model, $point)
println(
    "model, ", length(chains), " chains: ",
    round(median(trial).time / 1.0e3; digits = 2), " μs/point, ", trial.memory, " bytes/point",
)
trial = @benchmark sum(p -> unpolarized_intensity($model, p), $points)
println(length(points), " points: ", round(median(trial).time / 1.0e6; digits = 2), " ms")
