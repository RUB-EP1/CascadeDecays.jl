using CascadeDecays
using FourVectors
using HadronicLineshapes
using StaticArrays
using Test
using ThreeBodyDecays: NoRecoupling, RecouplingLS, VertexFunction

struct TestVertex <: AbstractVertex
    name::Symbol
end

@testset "DecayTopology" begin
    relation = [
        0 1
        0 1
        1 0
        1 -1
        -1 0
    ]
    topology = DecayTopology(relation; root = 5, finals = (1, 2, 3))

    @test validate_topology(topology)
    @test nlines(topology) == 5
    @test nvertices(topology) == 2
    @test nfinal(topology) == 3
    @test rootline(topology) == 5
    @test finallines(topology) == SVector(1, 2, 3)
    @test internal_lines(topology) == [4]
    @test has_canonical_line_order(topology)

    @test outgoing_lines(topology, 1) == [3, 4]
    @test child_lines(topology, 1) == [4, 3]
    @test final_descendants(topology, 4) == [1, 2]
    @test outgoing_lines(topology, 2) == [1, 2]
    @test child_lines(topology, 2) == [1, 2]

    @test produced_by(topology, 4) == 1
    @test consumed_by(topology, 4) == 2
    @test produced_by(topology, 5) === nothing
    @test consumed_by(topology, 3) === nothing

    @test bracket(topology) == "((1,2),3)"
end

@testset "CascadeSystem and kinematic routing" begin
    relation = [
        0 1
        0 1
        1 0
        1 -1
        -1 0
    ]
    topology = DecayTopology(relation; root = 5, finals = (1, 2, 3))
    system = CascadeSystem((0, 0, 0, 2), (0.1, 0.2, 0.3, 9.0))
    chain = DecayChain(
        topology,
        (ConstantLineshape(1.0),),
        (TestVertex(:mother), TestVertex(:isobar)),
        (4,),
        (2,),
    )
    x = CascadeKinematics(
        topology,
        system;
        internal_masses2 = (4.0,),
        vertex_angles = ((cosθ = 0.5, ϕ = 0.1), (cosθ = -0.2, ϕ = 0.3)),
    )

    @test line_masses2(topology, system, (4.0,)) == SVector(0.1, 0.2, 0.3, 4.0, 9.0)
    @test line_invariant(x, 4) == 4.0
    @test line_invariant(topology, x, (1, 2)) == 4.0

    @test vertex_masses2(topology, x, 1) == (9.0, 4.0, 0.3)
    @test vertex_masses2(topology, x, 2) == (4.0, 0.1, 0.2)
    @test vertex_masses2(topology, x, (1, 2)) == (4.0, 0.1, 0.2)
    @test line_two_js(chain, system) == SVector(0, 0, 0, 2, 2)
    @test vertex_spins(chain, system, 1) == (2, 2, 0)
    @test vertex_spins(chain, system, 2) == (2, 0, 0)
    @test vertex_helicities(topology, (0, 0, 0, 1, -1), 1) == (-1, 1, 0)
    @test vertex_helicities(topology, (0, 0, 0, 1, -1), 2) == (1, 0, 0)
    @test vertex_angles(x, 1) == (cosθ = 0.5, ϕ = 0.1)
    @test vertex_angles(topology, x, (1, 2)) == (cosθ = -0.2, ϕ = 0.3)

    @test_throws ArgumentError CascadeKinematics(
        topology,
        system;
        internal_masses2 = (),
        vertex_angles = ((cosθ = 0.5, ϕ = 0.1), (cosθ = -0.2, ϕ = 0.3)),
    )
end

@testset "DecayTopology validation" begin
    invalid_vertex = [
        0 1
        0 0
        1 0
        1 -1
        -1 0
    ]
    @test_throws ArgumentError DecayTopology(invalid_vertex; root = 5, finals = (1, 2, 3))

    invalid_final = [
        0 1
        0 1
        1 0
        1 -1
        -1 0
    ]
    @test_throws ArgumentError DecayTopology(invalid_final; root = 5, finals = (1, 2, 4))

    invalid_root = [
        0 1
        0 1
        1 0
        1 -1
        -1 0
    ]
    @test_throws ArgumentError DecayTopology(invalid_root; root = 1, finals = (1, 2, 3))
    @test_throws ArgumentError incoming_lines(DecayTopology(invalid_root; root = 5, finals = (1, 2, 3)), 3)
    @test_throws ArgumentError produced_by(DecayTopology(invalid_root; root = 5, finals = (1, 2, 3)), 6)
end

@testset "DecayChain payload mapping" begin
    relation = [
        0 1
        0 1
        1 0
        1 -1
        -1 0
    ]
    topology = DecayTopology(relation; root = 5, finals = (1, 2, 3))
    chain = DecayChain(
        topology,
        (ConstantLineshape(1.0),),
        (TestVertex(:mother), TestVertex(:isobar)),
        (4,),
        (2,),
    )

    @test chain.topology === topology
    @test nfinal(chain) == 3
    @test chain.propagator_two_js == SVector(2)
    @test bracket(chain) == "((1,2),3)"
    @test final_descendants(chain, 4) == [1, 2]
    @test chain.propagators[1](2.0) == 1.0

    @test_throws ArgumentError DecayChain(
        topology,
        (ConstantLineshape(1.0),),
        (TestVertex(:mother), TestVertex(:isobar)),
        (1,),
        (2,),
    )
    @test_throws ArgumentError DecayChain(
        topology,
        (ConstantLineshape(1.0), ConstantLineshape(2.0)),
        (TestVertex(:mother), TestVertex(:isobar)),
        (4,),
        (2, 4),
    )
end

@testset "Bracket-addressed DecayChain constructor" begin
    topology = DecayTopology((((1, 2), 3), 4))
    chain = DecayChain(
        topology;
        propagators = (
            ((1, 2), 3) => (two_j = 4, lineshape = ConstantLineshape(:R123)),
            (1, 2) => (two_j = 2, lineshape = ConstantLineshape(:R12)),
        ),
        vertices = (
            (1, 2) => TestVertex(:inner),
            (((1, 2), 3), 4) => TestVertex(:root),
            ((1, 2), 3) => TestVertex(:middle),
        ),
    )
    system = CascadeSystem((0, 0, 0, 0, 2), (0.1, 0.2, 0.3, 0.4, 9.0))

    @test chain.propagator_two_js == SVector(4, 2)
    @test chain.propagators[1](1.0) == :R123
    @test chain.propagators[2](1.0) == :R12
    @test chain.vertices[1].name == :root
    @test chain.vertices[2].name == :middle
    @test chain.vertices[3].name == :inner
    @test line_two_js(chain, system) == SVector(0, 0, 0, 0, 2, 4, 2)
    @test line_values(
        topology;
        finals = (0, 0, 0, 0),
        internals = ((1, 2) => 2, ((1, 2), 3) => -2),
        root = 0,
    ) == SVector(0, 0, 0, 0, 2, -2, 0)

    @test_throws ArgumentError DecayChain(
        topology;
        propagators = ((1, 2) => (two_j = 2, lineshape = ConstantLineshape(:R12)),),
        vertices = (
            (((1, 2), 3), 4) => TestVertex(:root),
            ((1, 2), 3) => TestVertex(:middle),
            (1, 2) => TestVertex(:inner),
        ),
    )
    @test_throws ArgumentError DecayChain(
        topology;
        propagators = (
            ((1, 2), 3) => (two_j = 4, lineshape = ConstantLineshape(:R123)),
            (1, 2) => (lineshape = ConstantLineshape(:R12),),
        ),
        vertices = (
            (((1, 2), 3), 4) => TestVertex(:root),
            ((1, 2), 3) => TestVertex(:middle),
            (1, 2) => TestVertex(:inner),
        ),
    )
end

@testset "Canonical routing on larger tree" begin
    relation = [
        0 0 1
        0 0 1
        0 1 0
        1 0 0
        0 1 -1
        1 -1 0
        -1 0 0
    ]
    topology = DecayTopology(relation; root = 7, finals = (1, 2, 3, 4))
    system = CascadeSystem((0, 0, 0, 0, 2), (0.1, 0.2, 0.3, 0.4, 9.0))
    x = CascadeKinematics(
        topology,
        system;
        internal_masses2 = (1.2, 2.3),
        vertex_angles = (
            (cosθ = 0.1, ϕ = 0.2),
            (cosθ = 0.3, ϕ = 0.4),
            (cosθ = 0.5, ϕ = 0.6),
        ),
    )
    chain = DecayChain(
        topology,
        (ConstantLineshape(11), ConstantLineshape(13)),
        (TestVertex(:root), TestVertex(:middle), TestVertex(:inner)),
        (6, 5),
        (4, 2),
    )

    @test has_canonical_line_order(topology)
    @test outgoing_lines(topology, 1) == [4, 6]
    @test child_lines(topology, 1) == [6, 4]
    @test bracket(topology) == "(((1,2),3),4)"

    @test vertex_masses2(topology, x, 1) == (9.0, 2.3, 0.4)
    @test vertex_masses2(topology, x, 2) == (2.3, 1.2, 0.3)
    @test vertex_masses2(topology, x, 3) == (1.2, 0.1, 0.2)
    @test vertex_masses2(topology, x, ((1, 2), 3)) == (2.3, 1.2, 0.3)
    @test line_two_js(chain, system) == SVector(0, 0, 0, 0, 2, 4, 2)
    @test vertex_spins(chain, system, 2) == (4, 2, 0)
    @test vertex_helicities(topology, (0, 0, 0, 0, 1, -2, 2), 2) == (-2, 1, 0)
    @test routed_propagator_product(chain, x) == 143
    @test line_invariant(topology, x, ((1, 2), 3)) == 2.3
    @test vertex_angles(topology, x, (((1, 2), 3), 4)) == (cosθ = 0.1, ϕ = 0.2)

    @test_throws ArgumentError line_invariant(x, 8)
    @test_throws ArgumentError vertex_angles(x, 4)
    @test_throws ArgumentError CascadeKinematics(
        topology,
        system;
        internal_masses2 = (1.2, 2.3),
        vertex_angles = ((θ = 0.1, ϕ = 0.2), (cosθ = 0.3, ϕ = 0.4), (cosθ = 0.5, ϕ = 0.6)),
    )
end

@testset "Bracket topology constructor" begin
    topology = DecayTopology((((1, 2), 3), 4))
    relation = [
        0 0 1
        0 0 1
        0 1 0
        1 0 0
        0 1 -1
        1 -1 0
        -1 0 0
    ]

    @test CascadeDecays.relation(topology) == SMatrix{7,3,Int,21}(relation)
    @test nfinal(topology) == 4
    @test nvertices(topology) == 3
    @test bracket(topology) == "(((1,2),3),4)"

    @test_throws ArgumentError DecayTopology(((1, 3), 4))
    @test_throws ArgumentError DecayTopology((1, (2, (3, 3))))
    @test_throws ArgumentError DecayTopology((1, 2, 3))
end

@testset "Topology-generated kinematics" begin
    topology = DecayTopology((((1, 2), 3), 4))

    pDminus = FourVector(0.8634762475578601, -0.2273640501540901, 0.5897962254778486; E = 2.1542368373711818)
    pD0 = FourVector(-0.41098561980779524, 0.4602903155362023, 0.1331926788204417; E = 1.9687832721903593)
    pKplus = FourVector(-0.4483316412801869, -0.23415599517164887, -0.7305348220308255; E = 1.0164784292725653)
    piplus = FourVector(-0.004158986333048991, 0.0012297296374225927, 0.00754591768614862; E = 0.13984149613322175)
    objs = (pD0, piplus, pDminus, pKplus)

    system = CascadeSystem((0, 0, 0, 0, 0), (mass.(objs) .^ 2..., mass(sum(objs))^2))
    programs = helicity_angle_programs(topology)
    x = cascade_kinematics(topology, system, objs)

    @test length(programs) == 3
    @test length.(programs) == (2, 3, 4)
    @test line_invariant(topology, x, (1, 2)) ≈ mass(pD0 + piplus)^2
    @test line_invariant(topology, x, ((1, 2), 3)) ≈ mass(pD0 + piplus + pDminus)^2
    @test line_invariant(topology, x, (((1, 2), 3), 4)) ≈ mass(sum(objs))^2
    @test isapprox(vertex_angles(topology, x, (1, 2)).cosθ, 0.8746538492596707; atol = 2e-10)
    @test isapprox(vertex_angles(topology, x, (1, 2)).ϕ, -2.84901364039537; atol = 2e-10)

    chain = DecayChain(
        topology;
        propagators = (
            (1, 2) => (two_j = 2, lineshape = ConstantLineshape(1.0 + 0.0im)),
            ((1, 2), 3) => (two_j = 2, lineshape = BreitWigner(4.039, 0.08)),
        ),
        vertices = (
            (((1, 2), 3), 4) => VertexFunction(RecouplingLS((2, 2))),
            ((1, 2), 3) => VertexFunction(RecouplingLS((0, 2))),
            (1, 2) => VertexFunction(RecouplingLS((2, 0))),
        ),
    )
    A = amplitude(chain, system, x, (0, 0, 0, 0, 0))
    @test isfinite(real(A))
    @test isfinite(imag(A))
end

@testset "Particle-2 helicity phase" begin
    masses2 = (1.0, 0.25, 0.25)
    angles = (cosθ = 1.0, ϕ = 0.0)
    spins = (0, 1, 1)

    positive_phase = VertexFunction(NoRecoupling(1, 1))
    negative_phase = VertexFunction(NoRecoupling(-1, -1))

    @test routed_vertex_amplitude(positive_phase, masses2, (0, 1, 1), spins, angles) == 1
    @test routed_vertex_amplitude(negative_phase, masses2, (0, -1, -1), spins, angles) == -1
end

include("threebody_compat_tests.jl")
