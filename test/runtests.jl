using CascadeDecays
using StaticArrays
using Test

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

    @test incoming_line(topology, 1) == 5
    @test outgoing_lines(topology, 1) == [3, 4]
    @test child_lines(topology, 1) == [4, 3]
    @test vertex_lines(topology, 1) == (5, 4, 3)
    @test final_descendants(topology, 4) == [1, 2]
    @test incoming_line(topology, 2) == 4
    @test outgoing_lines(topology, 2) == [1, 2]
    @test child_lines(topology, 2) == [1, 2]
    @test vertex_lines(topology, 2) == (4, 1, 2)

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
    system = CascadeSystem(
        (0, 0, 0, 2, 2),
        (0.1, 0.2, 0.3);
        root_mass2 = 9.0,
    )
    x = CascadeKinematics(
        topology,
        system;
        internal_masses2 = (4.0,),
        vertex_angles = ((cosθ = 0.5, ϕ = 0.1), (cosθ = -0.2, ϕ = 0.3)),
    )

    @test line_masses2(topology, system, (4.0,)) == SVector(0.1, 0.2, 0.3, 4.0, 9.0)
    @test x.line_masses2 == SVector(0.1, 0.2, 0.3, 4.0, 9.0)
    @test line_invariant(x, 4) == 4.0

    @test vertex_masses2(topology, x, 1) == (9.0, 4.0, 0.3)
    @test vertex_masses2(topology, x, 2) == (4.0, 0.1, 0.2)
    @test vertex_spins(topology, system, 1) == (2, 2, 0)
    @test vertex_spins(topology, system, 2) == (2, 0, 0)
    @test vertex_helicities(topology, (0, 0, 0, 1, -1), 1) == (-1, 1, 0)
    @test vertex_helicities(topology, (0, 0, 0, 1, -1), 2) == (1, 0, 0)
    @test vertex_angles(x, 1) == (cosθ = 0.5, ϕ = 0.1)

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
    )

    @test chain.topology === topology
    @test nfinal(chain) == 3
    @test propagating_lines(chain) == SVector(4)
    @test bracket(chain) == "((1,2),3)"
    @test vertex_lines(chain, 1) == (5, 4, 3)
    @test final_descendants(chain, 4) == [1, 2]
    @test chain.propagators[1](2.0) == 1.0

    @test_throws ArgumentError DecayChain(
        topology,
        (ConstantLineshape(1.0),),
        (TestVertex(:mother), TestVertex(:isobar)),
        (1,),
    )
    @test_throws ArgumentError DecayChain(
        topology,
        (ConstantLineshape(1.0), ConstantLineshape(2.0)),
        (TestVertex(:mother), TestVertex(:isobar)),
        (4,),
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
    system = CascadeSystem(
        (0, 0, 0, 0, 2, 4, 2),
        (0.1, 0.2, 0.3, 0.4);
        root_mass2 = 9.0,
    )
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
        (ConstantLineshape(:R123), ConstantLineshape(:R12)),
        (TestVertex(:root), TestVertex(:middle), TestVertex(:inner)),
        (6, 5),
    )

    @test has_canonical_line_order(topology)
    @test outgoing_lines(topology, 1) == [4, 6]
    @test child_lines(topology, 1) == [6, 4]
    @test vertex_lines(topology, 1) == (7, 6, 4)
    @test vertex_lines(topology, 2) == (6, 5, 3)
    @test vertex_lines(topology, 3) == (5, 1, 2)
    @test bracket(topology) == "(((1,2),3),4)"

    @test x.line_masses2 == SVector(0.1, 0.2, 0.3, 0.4, 1.2, 2.3, 9.0)
    @test vertex_masses2(topology, x, 1) == (9.0, 2.3, 0.4)
    @test vertex_masses2(topology, x, 2) == (2.3, 1.2, 0.3)
    @test vertex_masses2(topology, x, 3) == (1.2, 0.1, 0.2)
    @test vertex_spins(topology, system, 2) == (4, 2, 0)
    @test vertex_helicities(topology, (0, 0, 0, 0, 1, -2, 2), 2) == (-2, 1, 0)
    @test propagating_lines(chain) == SVector(6, 5)
    @test chain.propagators[1](line_invariant(x, propagating_lines(chain)[1])) == :R123
    @test chain.propagators[2](line_invariant(x, propagating_lines(chain)[2])) == :R12

    @test_throws ArgumentError line_invariant(x, 8)
    @test_throws ArgumentError vertex_angles(x, 4)
    @test_throws ArgumentError CascadeKinematics(
        topology,
        system;
        internal_masses2 = (1.2, 2.3),
        vertex_angles = ((θ = 0.1, ϕ = 0.2), (cosθ = 0.3, ϕ = 0.4), (cosθ = 0.5, ϕ = 0.6)),
    )
end
