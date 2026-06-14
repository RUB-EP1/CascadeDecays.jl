include(joinpath(@__DIR__, "..", "compat", "ThreeBodyCompat.jl"))
using .ThreeBodyCompat

using FourVectors
using ThreeBodyDecays:
    ThreeBodyMasses,
    ThreeBodyParities,
    ThreeBodySpins,
    ThreeBodySystem,
    aligned_four_vectors,
    cosθij,
    x2σs,
    @jp_str
import ThreeBodyDecays: DecayChainLS, amplitude as tbd_amplitude
import CascadeDecays:
    DecayChain,
    DecayTopology,
    CascadeDecay,
    CascadeKinematics,
    CascadeSystem,
    ConstantLineshape,
    CurrentFrame,
    HelicityRootFrame,
    KinematicTask,
    KinematicPoint,
    Propagator,
    RecouplingLS,
    SystemMasses,
    SystemSpins,
    Vertex,
    amplitude as cascade_amplitude,
    bracket,
    external_wigner_angles,
    kinematics_at,
    line_invariant,
    minimal_ls_decay_chain,
    vertex_angles,
    vertex_masses2

function _fourvector_from_tuple(p)
    return FourVector(Float64(p[1]), Float64(p[2]), Float64(p[3]); E = Float64(p[4]))
end

function _threebody_reference_point()
    ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0 = 7.7)
    σs = x2σs([0.5, 0.3], ms; k = 3)
    topology = three_body_topology(3)
    system = CascadeSystem(SystemSpins(0, 0, 0; two_h0 = 0), SystemMasses(ms))
    return (; ms, σs, topology, system)
end

@testset "three_body_topology resonance–spectator convention" begin
    @test bracket(three_body_topology(3)) == "((1,2),3)"
    @test bracket(three_body_topology(2)) == "((3,1),2)"
    @test bracket(DecayTopology(((1, 3), 2))) == "((1,3),2)"
    @test bracket(DecayTopology(((3, 1), 2))) == "((3,1),2)"
    @test DecayTopology(((1, 3), 2)) != DecayTopology(((3, 1), 2))
    @test three_body_k(three_body_topology(1)) == 1
    @test three_body_isobar_pair(1) == (2, 3)
    @test three_body_isobar_pair(2) == (3, 1)
    @test three_body_root_vertex(1) == ((2, 3), 1)
    @test three_body_isobar_pair(3) == (1, 2)
    @test three_body_root_vertex(3) == ((1, 2), 3)
end

@testset "ThreeBodyDecays aligned kinematics compatibility" begin
    (; ms, σs, topology, system) = _threebody_reference_point()
    expected_inner_cosθ = cosθij(σs, ms^2; k = 3)

    hand_x = CascadeKinematics(
        topology,
        system;
        internal_masses2 = (σs[3],),
        vertex_angles = ((cosθ = 1.0, ϕ = 0.0), (cosθ = expected_inner_cosθ, ϕ = 0.0)),
    )

    @test line_invariant(topology, hand_x, (1, 2)) ≈ σs[3]
    @test vertex_masses2(topology, hand_x, ((1, 2), 3)) == (ms.m0^2, σs[3], ms[3]^2)
    @test vertex_masses2(topology, hand_x, (1, 2)) == (σs[3], ms[1]^2, ms[2]^2)
    @test vertex_angles(topology, hand_x, ((1, 2), 3)) == (cosθ = 1.0, ϕ = 0.0)
    @test vertex_angles(topology, hand_x, (1, 2)) == (cosθ = expected_inner_cosθ, ϕ = 0.0)

    objs = _fourvector_from_tuple.(aligned_four_vectors(σs, ms; k = 3))
    generated_x = CascadeKinematics(topology, objs)

    @test line_invariant(topology, generated_x, (1, 2)) ≈ σs[3]
    @test vertex_angles(topology, generated_x, ((1, 2), 3)).cosθ ≈ 1.0
    @test vertex_angles(topology, generated_x, ((1, 2), 3)).ϕ ≈ 0.0
    @test vertex_angles(topology, generated_x, (1, 2)).cosθ ≈ expected_inner_cosθ
    @test vertex_angles(topology, generated_x, (1, 2)).ϕ ≈ 0.0
    @test three_body_refζs(topology) == (3, 3, 3, 3)
end

@testset "Rotated aligned root kinematics diagnostic" begin
    (; ms, σs, topology, system) = _threebody_reference_point()
    objs = _fourvector_from_tuple.(aligned_four_vectors(σs, ms; k = 3))

    θ = deg2rad(10.0)
    ϕ = deg2rad(20.0)
    rotated_objs = map(p -> p |> Ry(θ) |> Rz(ϕ), objs)
    generated_x = CascadeKinematics(topology, rotated_objs; initial_frame = CurrentFrame())
    root_angles = vertex_angles(topology, generated_x, ((1, 2), 3))

    @test root_angles.cosθ ≈ cos(θ)
    @test root_angles.ϕ ≈ ϕ
end

@testset "CurrentFrame kinematics preserve aligned rotation convention" begin
    (; ms, σs, topology, system) = _threebody_reference_point()
    objs_cm = _fourvector_from_tuple.(aligned_four_vectors(σs, ms; k = 3))
    event_transform(p) = p |> Rz(0.5) |> Ry(0.3) |> Rz(0.4)
    objs = Tuple(event_transform(p) for p in objs_cm)

    generated_x = CascadeKinematics(topology, objs; initial_frame = CurrentFrame())
    root_angles = vertex_angles(topology, generated_x, ((1, 2), 3))
    inner_angles = vertex_angles(topology, generated_x, (1, 2))
    po = reference_plane_orientation(topology, objs; initial_frame = CurrentFrame())
    po_from_root_frame = reference_plane_orientation(topology, objs; initial_frame = HelicityRootFrame())

    @test root_angles.ϕ ≈ 0.4
    @test root_angles.cosθ ≈ cos(0.3)
    @test inner_angles.ϕ ≈ 0.5
    @test inner_angles.cosθ ≈ cosθij(σs, ms^2; k = 3)
    @test po.α ≈ 0.4
    @test po.cosβ ≈ cos(0.3)
    @test po.γ ≈ 0.5
    @test po_from_root_frame.α ≈ po.α
    @test po_from_root_frame.cosβ ≈ po.cosβ
    @test po_from_root_frame.γ ≈ po.γ
end

@testset "HelicityRootFrame kinematics are invariant to later boosts" begin
    (; ms, σs, topology, system) = _threebody_reference_point()
    objs_cm = _fourvector_from_tuple.(aligned_four_vectors(σs, ms; k = 3))
    event_transform(p) = p |> Rz(0.5) |> Ry(0.3) |> Rz(0.4) |> Bz(3.3) |> Ry(0.1) |> Rz(0.1)
    objs = Tuple(event_transform(p) for p in objs_cm)

    generated_x = CascadeKinematics(topology, objs; initial_frame = HelicityRootFrame())
    root_angles = vertex_angles(topology, generated_x, ((1, 2), 3))
    inner_angles = vertex_angles(topology, generated_x, (1, 2))

    @test root_angles.ϕ ≈ 0.4
    @test root_angles.cosθ ≈ cos(0.3)
    @test inner_angles.ϕ ≈ 0.5
    @test inner_angles.cosθ ≈ cosθij(σs, ms^2; k = 3)
end

@testset "HelicityRootFrame ignores numerical rest-frame residuals" begin
    (; ms, σs, topology, system) = _threebody_reference_point()
    objs_cm = _fourvector_from_tuple.(aligned_four_vectors(σs, ms; k = 3))
    rotated_objs = Tuple(p |> Ry(0.6) |> Rz(0.5) for p in objs_cm)
    noisy_objs = (
        FourVector(
            rotated_objs[1].px + 1e-13,
            rotated_objs[1].py + 2e-13,
            rotated_objs[1].pz - 1e-13;
            E = rotated_objs[1].E,
        ),
        rotated_objs[2],
        rotated_objs[3],
    )

    generated_x = CascadeKinematics(topology, noisy_objs; initial_frame = HelicityRootFrame())
    root_angles = vertex_angles(topology, generated_x, ((1, 2), 3))
    inner_angles = vertex_angles(topology, generated_x, (1, 2))

    @test root_angles.ϕ ≈ 0.5
    @test root_angles.cosθ ≈ cos(0.6)
    @test inner_angles.ϕ ≈ 0.0 atol = 1e-12
    @test inner_angles.cosθ ≈ cosθij(σs, ms^2; k = 3)
end

function _threebody_cascade_chain(;
    jp = "0+",
    two_js = ThreeBodySpins(0, 0, 0; two_h0 = 0),
    Xlineshape = σ -> 1.0 + 0.0im,
)
    ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0 = 7.7)
    σs = x2σs([0.5, 0.3], ms; k = 3)
    Ps = ThreeBodyParities('+', '+', '+'; P0 = '+')
    dc = DecayChainLS(; k = 3, Xlineshape, jp, Ps, tbs = ThreeBodySystem(ms, two_js))

    chain_topology = three_body_topology(3)
    reference_topology = default_three_body_reference_topology()
    system = CascadeSystem(
        SystemSpins(two_js.two_h1, two_js.two_h2, two_js.two_h3; two_h0 = two_js.two_h0),
        SystemMasses(ms),
    )
    isobar = three_body_isobar_pair(3)
    root = three_body_root_vertex(3)
    chain = DecayChain(
        chain_topology;
        propagators = (isobar => Propagator(dc.two_j, Xlineshape),),
        vertices = (
            root => Vertex(RecouplingLS(dc.HRk.h.two_ls)),
            isobar => Vertex(RecouplingLS(dc.Hij.h.two_ls)),
        ),
    )
    objs = _fourvector_from_tuple.(aligned_four_vectors(σs, ms; k = 3))
    task = KinematicTask(
        (three_body_topology(1), chain_topology);
        reference_topology = reference_topology,
        wigner_finals = (1, 2, 3),
    )
    point = KinematicPoint(task, objs)
    x = kinematics_at(point, chain_topology)
    return (; ms, σs, dc, chain, reference_topology, system, objs, x, point)
end

@testset "ThreeBodyDecays helicity amplitude |A|² (scalar chain)" begin
    (; σs, dc, chain, reference_topology, system, objs, point) =
        _threebody_cascade_chain(; jp = "0+", two_js = ThreeBodySpins(0, 0, 0; two_h0 = 0))

    po = reference_plane_orientation(chain.topology, objs)
    refζs = three_body_refζs(chain.topology)
    A_tbd = tbd_amplitude(dc, po, σs; refζs = refζs)
    A_cascade = cascade_amplitude(chain, system, point)

    @test refζs == (3, 3, 3, 3)
    @test size(A_cascade) == size(A_tbd)
    @test isapprox(A_cascade, A_tbd)
    @test sum(abs2, A_cascade) ≈ sum(abs2, A_tbd)
    @test sum(abs2, A_cascade) ≈ 1.0
end

@testset "generated root and isobar vertex angles" begin
    (; chain, x) = _threebody_cascade_chain()
    root = vertex_angles(chain.topology, x, three_body_root_vertex(3))
    inner = vertex_angles(chain.topology, x, three_body_isobar_pair(3))
    @test root.ϕ ≈ 0.0
    @test root.cosθ ≈ 1.0
    @test inner.ϕ ≈ 0.0
end

@testset "Identity reference topology gives trivial relative Wigner angles" begin
    ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0 = 7.7)
    σs = x2σs([0.5, 0.3], ms; k = 3)
    topology = three_body_topology(3)
    objs = _fourvector_from_tuple.(aligned_four_vectors(σs, ms; k = 3))
    angles = external_wigner_angles(topology, topology, objs)
    for ang in angles
        @test ang == (ϕ = 0.0, θ = 0.0, ψ = 0.0)
    end
end

@testset "CascadeDecay container amplitude routing" begin
    (; chain, reference_topology, system, point) = _threebody_cascade_chain()
    cascade = CascadeDecay((chain,), system, reference_topology; couplings = (1.0,))
    @test cascade_amplitude(cascade, point) ≈ cascade_amplitude(chain, system, point)
end
