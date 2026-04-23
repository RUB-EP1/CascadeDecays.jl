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

    @test incoming_line(topology, 1) == 5
    @test outgoing_lines(topology, 1) == [3, 4]
    @test incoming_line(topology, 2) == 4
    @test outgoing_lines(topology, 2) == [1, 2]

    @test produced_by(topology, 4) == 1
    @test consumed_by(topology, 4) == 2
    @test produced_by(topology, 5) === nothing
    @test consumed_by(topology, 3) === nothing

    @test bracket(topology) == "((1,2),3)"
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
    @test chain.propagators[1](2.0) == 1.0

    @test_throws ArgumentError DecayChain(
        topology,
        (ConstantLineshape(1.0),),
        (TestVertex(:mother), TestVertex(:isobar)),
        (1,),
    )
end
