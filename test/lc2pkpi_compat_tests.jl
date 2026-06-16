# Λc⁺ → p K⁻ π⁺ single-chain cross-check against ThreeBodyDecays.jl.
include(joinpath(@__DIR__, "..", "compat", "ThreeBodyCompat.jl"))
using .ThreeBodyCompat

using FourVectors
using ThreeBodyDecays:
    ThreeBodyMasses,
    ThreeBodyParities,
    ThreeBodySpins,
    ThreeBodySystem,
    aligned_four_vectors,
    Invariants,
    x2σs,
    @jp_str
import ThreeBodyDecays: DecayChainLS, amplitude as tbd_amplitude
import CascadeDecays:
    CascadeSystem,
    ConstantLineshape,
    KinematicTask,
    KinematicPoint,
    Propagator,
    SystemMasses,
    SystemSpinParities,
    amplitude as cascade_amplitude,
    minimal_ls_decay_chain

const LC_UNIT_LINESHAPE = ConstantLineshape(1.0 + 0.0im)
const LC_MASSES = (
    m0 = 2.28646,
    p = 0.93827208816,
    K = 0.493677,
    π = 0.13957039,
)

function _lc_fourvector_from_tuple(p)
    return FourVector(Float64(p[1]), Float64(p[2]), Float64(p[3]); E = Float64(p[4]))
end

@testset "Lc2pKpi Lambda chain matches ThreeBodyDecays in top3 basis" begin
    ms = ThreeBodyMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0)
    tbs = ThreeBodySystem(ms, ThreeBodySpins(1, 0, 0; two_h0 = 1))
    Ps = ThreeBodyParities('+', '-', '-'; P0 = '+')
    unit_X(σ) = LC_UNIT_LINESHAPE(σ)

    tbd_chain = DecayChainLS(; k = 3, Xlineshape = unit_X, jp = jp"3/2-", Ps, tbs)
    system = CascadeSystem(
        SystemSpinParities("1/2+", "0-", "0-"; jp0 = "1/2+"),
        SystemMasses(LC_MASSES.p, LC_MASSES.K, LC_MASSES.π; m0 = LC_MASSES.m0),
    )
    chain = minimal_ls_decay_chain(
        three_body_topology(3),
        system,
        ((three_body_isobar_pair(3) => Propagator(jp"3/2-", LC_UNIT_LINESHAPE)),),
    )

    σs = x2σs([0.42, 0.31], ms; k = 3)
    objs = _lc_fourvector_from_tuple.(aligned_four_vectors(σs, ms; k = 3))
    task = KinematicTask((chain.topology,); reference_topology = chain.topology)
    point = KinematicPoint(task, objs)
    σtbd = Invariants(σs.σ1, σs.σ2, σs.σ3)
    po = reference_plane_orientation(chain.topology, objs)

    A_cascade = cascade_amplitude(chain, system, point)
    A_tbd = tbd_amplitude(tbd_chain, po, σtbd; refζs = three_body_refζs(chain.topology))

    @test size(A_cascade) == (2, 1, 1, 2)
    @test A_cascade ≈ A_tbd
    @test sum(abs2, A_cascade) ≈ sum(abs2, A_tbd)
end
