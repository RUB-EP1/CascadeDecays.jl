using CascadeDecays
import CascadeDecays: possible_ls_more, UndefinedParity
using FourVectors
using HadronicLineshapes
using StaticArrays
using Test
using ThreeBodyDecays:
    @jp_str,
    NoRecoupling,
    RecouplingLS,
    SpinParity,
    ThreeBodyMasses,
    ThreeBodySpins,
    aligned_four_vectors,
    cosθij,
    x2σs

struct TestVertex <: AbstractVertex
    name::Symbol
end

@testset "Kinematic quantum-number types" begin
    ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0 = 7.7)
    @test SystemMasses(ms) == SystemMasses(1.1, 2.2, 3.3; m0 = 7.7)
    @test SystemSpins(ThreeBodySpins(0, 0, 0; two_h0 = 0)) == SystemSpins(0, 0, 0; two_h0 = 0)

    system = CascadeSystem(SystemSpins(0, 0, 0; two_h0 = 2), SystemMasses(1, 2, 3; m0 = 3))
    parity_system = CascadeSystem(SystemSpinParities(system.quantum, '+', '+', '+'; P0 = '+'), system.masses)
    @test parity_system.quantum.parities == SystemParities('+', '+', '+'; P0 = '+')

    from_strings = SystemSpinParities("1+", "0-", "0-"; jp0 = "1±")
    from_jp = SystemSpinParities(jp"1+", jp"0-", jp"0-"; jp0 = jp"1±")
    @test from_strings.spins == SystemSpins(2, 0, 0; two_h0 = 2)
    @test from_strings.parities == SystemParities('+', '-', '-'; P0 = UndefinedParity)
    @test from_jp == from_strings
    @test final_two_js(parity_system) == final_two_js(system)
    @test root_two_j(parity_system) == root_two_j(system)
    @test root_two_j(system) == 2
    @test root_mass(system) == 3

    @test_throws ArgumentError SystemSpins(0, 0)
    spins = SystemSpins(0, 0, 0, 0; two_h0 = 0)
    @test_throws ArgumentError SystemHelicities(spins, 0, 0, 0, 2; two_h0 = 0)
    @test SystemHelicities(spins, 0, 0, 0, 0; two_h0 = 0) ==
        SystemHelicities(0, 0, 0, 0; two_h0 = 0)
    @test SystemSpinProjections === SystemSpins
    @test_throws ArgumentError SystemHelicities(0, 0)
    @test_throws ArgumentError SystemHelicities(spins, 0, 0, 0)
    @test_throws ArgumentError SystemMasses(; m0 = 1.0)
end

@testset "DecayTopology" begin
    topology = DecayTopology(((1, 2), 3))

    @test nlines(topology) == 5
    @test nvertices(topology) == 2
    @test nfinal(topology) == 3
    @test root_line_ind(topology) == 5
    @test final_line_inds(topology) == SVector(1, 2, 3)
    @test internal_line_inds(topology) == [4]
    @test has_canonical_line_order(topology)

    @test outgoing_line_inds(topology, 1) == [3, 4]
    @test child_line_inds(topology, 1) == SVector(4, 3)
    @test final_descendants(topology, 4) == [1, 2]
    @test outgoing_line_inds(topology, 2) == [1, 2]
    @test child_line_inds(topology, 2) == [1, 2]

    @test produced_by(topology, 4) == 1
    @test consumed_by(topology, 4) == 2
    @test produced_by(topology, 5) === nothing
    @test consumed_by(topology, 3) === nothing

    @test bracket_notation(topology) == "((1,2),3)"
end

@testset "CascadeSystem and kinematic routing" begin
    topology = DecayTopology(((1, 2), 3))
    system = CascadeSystem(SystemSpins(0, 0, 0; two_h0 = 2), SystemMasses(1, 2, 3; m0 = 3))
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

    @test line_masses2(topology, system, (4.0,)) == SVector(1.0, 4.0, 9.0, 4.0, 9.0)
    @test line_invariant(x, 4) == 4.0
    @test line_invariant(topology, x, (1, 2)) == 4.0

    @test vertex_masses2(topology, x, 1) == (9.0, 4.0, 9.0)
    @test vertex_masses2(topology, x, 2) == (4.0, 1.0, 4.0)
    @test vertex_masses2(topology, x, (1, 2)) == (4.0, 1.0, 4.0)
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

function _fourvector_from_tuple(p)
    return FourVector(Float64(p[1]), Float64(p[2]), Float64(p[3]); E = Float64(p[4]))
end

@testset "KinematicTask from four-vectors" begin
    ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0 = 7.7)
    σs = x2σs([0.5, 0.3], ms; k = 3)
    ref_topology = DecayTopology(((1, 2), 3))
    alt_topology = DecayTopology(((3, 1), 2))

    objs_cm = Tuple(_fourvector_from_tuple(p) for p in aligned_four_vectors(σs, ms; k = 3))
    event_transform(p) = p |> Rz(0.5) |> Ry(0.3) |> Rz(0.4)
    objs = Tuple(event_transform(p) for p in objs_cm)

    task = KinematicTask(
        (ref_topology, alt_topology);
        reference_topology = ref_topology,
        wigner_finals = (1, 3),
        initial_frame = CurrentFrame(),
    )
    point = KinematicPoint(task, objs)
    x_ref = kinematics_at(point, ref_topology)

    @test length(task.programs) == 2
    @test length(task.programs[1].vertex_programs) == nvertices(ref_topology)
    program_show = sprint(show, task.programs[1].vertex_programs)
    @test occursin("VertexPrograms with 2 measurements", program_show)
    @test occursin("vertex ((1,2),3):", program_show)
    @test occursin("vertex (1,2):", program_show)
    @test bracket_notation(task.reference_topology) == "((1,2),3)"
    @test line_invariant(x_ref, 1) ≈ mass(objs[1])^2
    @test line_invariant(ref_topology, x_ref, (1, 2)) ≈ σs[3]
    @test line_invariant(x_ref, ref_topology, (1, 2)) ≈ σs[3]
    @test vertex_angles(x_ref, ref_topology, ((1, 2), 3)).ϕ ≈ 0.4
    @test vertex_angles(x_ref, ref_topology, ((1, 2), 3)).cosθ ≈ cos(0.3)
    @test vertex_angles(x_ref, ref_topology, (1, 2)).ϕ ≈ 0.5
    @test vertex_angles(x_ref, ref_topology, (1, 2)).cosθ ≈ cosθij(σs, ms^2; k = 3)

    ref_alignments = alignment_angles_at(point, ref_topology)
    alt_alignments = alignment_angles_at(point, alt_topology)
    @test length(ref_alignments) == 3
    @test length(alt_alignments) == 3
    @test ref_alignments[1] == (α = 0.0, cosβ = 1.0, γ = 0.0)
    @test ref_alignments[2] == (α = 0.0, cosβ = 1.0, γ = 0.0)
    @test ref_alignments[3] == (α = 0.0, cosβ = 1.0, γ = 0.0)
    @test alt_alignments[2] == (α = 0.0, cosβ = 1.0, γ = 0.0)
    @test alt_alignments[1] != (α = 0.0, cosβ = 1.0, γ = 0.0)
    @test alt_alignments[3] != (α = 0.0, cosβ = 1.0, γ = 0.0)

    boosted_objs = Tuple(p |> Bz(1.5) |> Ry(0.1) |> Rz(0.2) for p in objs)
    x_helicity = CascadeKinematics(ref_topology, boosted_objs; initial_frame = HelicityRootFrame())
    @test vertex_angles(x_helicity, ref_topology, ((1, 2), 3)).ϕ ≈ 0.4
    @test vertex_angles(x_helicity, ref_topology, ((1, 2), 3)).cosθ ≈ cos(0.3)
end

@testset "DecayChain payload mapping" begin
    topology = DecayTopology(((1, 2), 3))
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
    @test bracket_notation(chain) == "((1,2),3)"
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

@testset "CascadeDecay constructor" begin
    topology = DecayTopology(((1, 2), 3))
    system = CascadeSystem(SystemSpins(0, 0, 0; two_h0 = 0), SystemMasses(1, 2, 3; m0 = 3))
    chain = DecayChain(
        topology,
        (ConstantLineshape(1.0),),
        (TestVertex(:mother), TestVertex(:isobar)),
        (4,),
        (0,),
    )

    cascade = CascadeDecay((chain,), system, topology)
    @test cascade.chains === (chain,)
    @test cascade_system(cascade) === system
    @test reference_topology(cascade) === topology
    @test couplings(cascade) == (1.0 + 0.0im,)
    @test cascade.names == ("chain_1",)

    weighted = CascadeDecay((chain, chain), system, topology; couplings = (2, 3.0im))
    @test weighted.chains === (chain, chain)
    @test couplings(weighted) == (2.0 + 0.0im, 0.0 + 3.0im)

    amp_masses = ThreeBodyMasses(1.0, 1.0, 1.0; m0 = 3.5)
    amp_σs = x2σs([0.45, 0.35], amp_masses; k = 3)
    amp_objs = Tuple(_fourvector_from_tuple(p) for p in aligned_four_vectors(amp_σs, amp_masses; k = 3))
    amp_topology = DecayTopology(((1, 2), 3))
    amp_system = CascadeSystem(
        SystemSpins(0, 0, 0; two_h0 = 0),
        SystemMasses(amp_masses),
    )
    amp_chain = DecayChain(
        amp_topology;
        propagators = ((1, 2) => Propagator(0, ConstantLineshape(2.0)),),
        vertices = (
            (((1, 2), 3) => Vertex(NoRecoupling(0, 0))),
            ((1, 2) => Vertex(NoRecoupling(0, 0))),
        ),
    )
    amp_point = KinematicPoint(KinematicTask((amp_topology,); reference_topology = amp_topology), amp_objs)
    amp_cascade = CascadeDecay((amp_chain, amp_chain), amp_system, amp_topology; couplings = (2, 3.0im))
    @test amplitude(amp_cascade, amp_point) ≈ (2.0 + 3.0im) .* amplitude(amp_chain, amp_system, amp_point)
    @test unpolarized_intensity(amp_cascade, amp_point) ≈ sum(abs2, amplitude(amp_cascade, amp_point))

    @test_throws ArgumentError CascadeDecay((), system, topology)
    @test_throws ArgumentError CascadeDecay((chain,), system, topology; couplings = (1, 2))

    vector_cascade = CascadeDecay([chain], system, topology)
    @test vector_cascade.chains === (chain,)

    alt_topology = DecayTopology(((2, 3), 1))
    alt_chain = DecayChain(
        alt_topology,
        (ConstantLineshape(1.0),),
        (TestVertex(:mother), TestVertex(:isobar)),
        (4,),
        (0,),
    )
    mixed = CascadeDecay((chain, alt_chain), system, topology)
    @test mixed.chains[1].topology === topology
    @test mixed.chains[2].topology === alt_topology

    wrong_system = CascadeSystem(SystemSpins(0, 0; two_h0 = 0), SystemMasses(1, 2; m0 = 3))
    @test_throws ArgumentError CascadeDecay((chain,), wrong_system, topology)
end

@testset "CascadeDecay indexing, names, and merge" begin
    topology = DecayTopology(((1, 2), 3))
    system = CascadeSystem(SystemSpins(0, 0, 0; two_h0 = 0), SystemMasses(1, 2, 3; m0 = 3))
    chain_a = DecayChain(
        topology;
        propagators = ((1, 2) => Propagator(0, ConstantLineshape(1.0)),),
        vertices = (
            (((1, 2), 3) => Vertex(RecouplingLS((0, 0)))),
            ((1, 2) => Vertex(RecouplingLS((0, 0)))),
        ),
    )
    chain_b = DecayChain(
        topology;
        propagators = ((1, 2) => Propagator(0, ConstantLineshape(2.0)),),
        vertices = (
            (((1, 2), 3) => Vertex(RecouplingLS((2, 0)))),
            ((1, 2) => Vertex(RecouplingLS((0, 0)))),
        ),
    )

    model = CascadeDecay(
        (chain_a, chain_b),
        system,
        topology;
        couplings = (1.0, 2.0im),
        names = ("a1->rhopi", "a1->kk"),
    )
    @test model.names == ("a1->rhopi", "a1->kk")
    @test length(model) == 2
    @test collect(model)[1].name == "a1->rhopi"
    @test model[1].names == ("a1->rhopi",)
    @test model[1].couplings == (1.0,)
    @test model[[2]].names == ("a1->kk",)
    @test model[model.names .== "a1->rhopi"].names == ("a1->rhopi",)
    @test model["a1->kk"].chains === (chain_b,)

    default_names = CascadeDecay((chain_a, chain_b), system, topology).names
    @test default_names == ("chain_1", "chain_2")

    show_text = sprint(show, model)
    @test occursin("a1->rhopi", show_text)
    @test occursin("coupling", show_text)
    @test occursin("topology", show_text)
    @test occursin("((1,2),3)", show_text)

    part_a = model[model.names .== "a1->rhopi"]
    part_b = model[model.names .== "a1->kk"]
    combined = merge(part_a, part_b)
    @test combined.names == model.names
    @test combined.couplings == model.couplings
    @test combined.chains === model.chains

    amp_masses = ThreeBodyMasses(1.0, 1.0, 1.0; m0 = 3.5)
    amp_σs = x2σs([0.45, 0.35], amp_masses; k = 3)
    amp_objs = Tuple(_fourvector_from_tuple(p) for p in aligned_four_vectors(amp_σs, amp_masses; k = 3))
    amp_point = KinematicPoint(KinematicTask((topology,); reference_topology = topology), amp_objs)
    @test amplitude(model[[1]], amp_point) ≈ amplitude(part_a, amp_point)
    @test amplitude(model[[2]], amp_point) ≈ 2.0im * amplitude(chain_b, system, kinematics_at(amp_point, topology))

    @test_throws ArgumentError model[fill(false, length(model))]
    @test_throws KeyError model["missing"]
    @test_throws ArgumentError merge(part_a, part_a)
    @test_throws ArgumentError CascadeDecay(
        (chain_a, chain_b),
        system,
        topology;
        names = ("only-one",),
    )
end

@testset "Bracket-addressed DecayChain constructor" begin
    topology = DecayTopology((((1, 2), 3), 4))
    chain = DecayChain(
        topology;
        propagators = (
            ((1, 2), 3) => Propagator(4, ConstantLineshape(:R123)),
            (1, 2) => Propagator(2, ConstantLineshape(:R12)),
        ),
        vertices = (
            (1, 2) => TestVertex(:inner),
            (((1, 2), 3), 4) => TestVertex(:root),
            ((1, 2), 3) => TestVertex(:middle),
        ),
    )
    system = CascadeSystem(SystemSpins(0, 0, 0, 0; two_h0 = 2), SystemMasses(1, 2, 3, 4; m0 = 3))

    @test chain.propagator_two_js == SVector(4, 2)
    @test chain.propagators[1](1.0) == :R123
    @test chain.propagators[2](1.0) == :R12
    @test chain.vertices[1].name == :root
    @test chain.vertices[2].name == :middle
    @test chain.vertices[3].name == :inner
    @test line_two_js(chain, system) == SVector(0, 0, 0, 0, 2, 4, 2)

    jp_chain = DecayChain(
        topology;
        propagators = (
            ((1, 2), 3) => Propagator(SpinParity(4, '+'), ConstantLineshape(:R123)),
            (1, 2) => Propagator(SpinParity(2, '-'), ConstantLineshape(:R12)),
        ),
        vertices = (
            (((1, 2), 3), 4) => TestVertex(:root),
            ((1, 2), 3) => TestVertex(:middle),
            (1, 2) => TestVertex(:inner),
        ),
    )
    @test jp_chain.propagator_two_js == SVector(4, 2)
    @test line_values(
        topology;
        finals = (0, 0, 0, 0),
        internals = ((1, 2) => 2, ((1, 2), 3) => -2),
        root = 0,
    ) == SVector(0, 0, 0, 0, 2, -2, 0)

    @test_throws ArgumentError DecayChain(
        topology;
        propagators = ((1, 2) => Propagator(2, ConstantLineshape(:R12)),),
        vertices = (
            (((1, 2), 3), 4) => TestVertex(:root),
            ((1, 2), 3) => TestVertex(:middle),
            (1, 2) => TestVertex(:inner),
        ),
    )
    @test_throws TypeError DecayChain(
        topology;
        propagators = (
            ((1, 2), 3) => Propagator(4, ConstantLineshape(:R123)),
            (1, 2) => ConstantLineshape(:R12),
        ),
        vertices = (
            (((1, 2), 3), 4) => TestVertex(:root),
            ((1, 2), 3) => TestVertex(:middle),
            (1, 2) => TestVertex(:inner),
        ),
    )
end

@testset "Canonical routing on larger tree" begin
    topology = DecayTopology((((1, 2), 3), 4))
    system = CascadeSystem(SystemSpins(0, 0, 0, 0; two_h0 = 2), SystemMasses(1, 2, 3, 4; m0 = 3))
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
    @test outgoing_line_inds(topology, 1) == [4, 6]
    @test child_line_inds(topology, 1) == [6, 4]
    @test bracket_notation(topology) == "(((1,2),3),4)"

    @test vertex_masses2(topology, x, 1) == (9.0, 2.3, 16.0)
    @test vertex_masses2(topology, x, 2) == (2.3, 1.2, 9.0)
    @test vertex_masses2(topology, x, 3) == (1.2, 1.0, 4.0)
    @test vertex_masses2(topology, x, ((1, 2), 3)) == (2.3, 1.2, 9.0)
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

    @test CascadeDecays.relation(topology) == SMatrix{7, 3, Int, 21}(relation)
    @test nfinal(topology) == 4
    @test nvertices(topology) == 3
    @test bracket_notation(topology) == "(((1,2),3),4)"

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

    system = CascadeSystem(
        SystemSpins(0, 0, 0, 0; two_h0 = 0),
        SystemMasses(mass.(objs)...; m0 = mass(sum(objs))),
    )
    programs = CascadeDecays.helicity_angle_programs(topology)
    x = CascadeKinematics(topology, objs)

    @test length(programs) == 3
    @test length.(programs) == (2, 3, 4)
    @test line_invariant(topology, x, (1, 2)) ≈ mass(pD0 + piplus)^2
    @test line_invariant(topology, x, ((1, 2), 3)) ≈ mass(pD0 + piplus + pDminus)^2
    @test line_invariant(topology, x, (((1, 2), 3), 4)) ≈ mass(sum(objs))^2
    @test isapprox(vertex_angles(topology, x, (1, 2)).cosθ, 0.8746538492596707; atol = 2.0e-10)
    @test isapprox(vertex_angles(topology, x, (1, 2)).ϕ, -2.84901364039537; atol = 2.0e-10)

    chain = DecayChain(
        topology;
        propagators = (
            (1, 2) => Propagator(2, ConstantLineshape(1.0 + 0.0im)),
            ((1, 2), 3) => Propagator(2, BreitWigner(4.039, 0.08)),
        ),
        vertices = (
            (((1, 2), 3), 4) => Vertex(RecouplingLS((2, 2))),
            ((1, 2), 3) => Vertex(RecouplingLS((0, 2))),
            (1, 2) => Vertex(RecouplingLS((2, 0))),
        ),
    )
    external_two_λs = SystemHelicities(system.quantum, 0, 0, 0, 0; two_h0 = 0)
    A_full = amplitude(chain, system, x)
    @test size(A_full) == (1, 1, 1, 1, 1)
    A = amplitude(chain, system, x, external_two_λs)
    @test A == A_full[1, 1, 1, 1, 1]
    @test isfinite(real(A))
    @test isfinite(imag(A))
    @test_throws MethodError amplitude(chain, system, x, (0, 0, 0, 0, 0))
end

@testset "Particle-2 helicity phase" begin
    masses2 = (1.0, 0.25, 0.25)
    angles = (cosθ = 1.0, ϕ = 0.0)
    spins = (0, 1, 1)

    positive_phase = Vertex(NoRecoupling(1, 1))
    negative_phase = Vertex(NoRecoupling(-1, -1))

    @test routed_vertex_amplitude(positive_phase, masses2, (0, 1, 1), spins, angles) == 1
    @test routed_vertex_amplitude(negative_phase, masses2, (0, -1, -1), spins, angles) == -1
end

@testset "KinematicPoint amplitude dispatch" begin
    ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0 = 3.5)
    σs = x2σs([0.45, 0.35], ms; k = 3)
    objs = Tuple(_fourvector_from_tuple(p) for p in aligned_four_vectors(σs, ms; k = 3))
    topology = DecayTopology(((1, 2), 3))
    system = CascadeSystem(
        SystemSpins(0, 0, 0; two_h0 = 0),
        SystemMasses(ms),
    )
    chain = DecayChain(
        topology;
        propagators = ((1, 2) => Propagator(0, ConstantLineshape(2.0)),),
        vertices = (
            (((1, 2), 3) => Vertex(NoRecoupling(0, 0))),
            ((1, 2) => Vertex(NoRecoupling(0, 0))),
        ),
    )
    point = KinematicPoint(KinematicTask((topology,); reference_topology = topology), objs)
    x = kinematics_at(point, topology)

    @test amplitude(chain, system, point) == amplitude(chain, system, x)
end

@testset "LS decay-chain builders" begin
    topology = DecayTopology(((1, 2), 3))
    system = CascadeSystem(
        SystemSpins(0, 0, 0; two_h0 = 0),
        SystemMasses(1.0, 1.0, 1.0; m0 = 3.0),
    )
    weak_system = CascadeSystem(
        SystemSpinParities(system.quantum, '+', '+', '+'; P0 = '+'),
        system.masses,
    )
    propagators = (
        (1, 2) => Propagator(SpinParity(0, '+'), ConstantLineshape(1.0)),
    )

    @test minimal_vertex_couplings(topology, weak_system, propagators) == (
        (((1, 2), 3) => (0, 0)),
        ((1, 2) => (0, 0)),
    )
    @test_throws MethodError possible_vertex_couplings(
        topology,
        weak_system,
        ((1, 2) => Propagator(0, ConstantLineshape(1.0)),),
    )

    minimal = minimal_ls_decay_chain(topology, weak_system, propagators)
    @test minimal.propagator_two_js == SVector(0)
    @test all(v -> v.h isa RecouplingLS, minimal.vertices)

    all_chains = all_ls_decay_chains(topology, weak_system, propagators)
    @test length(all_chains) == 1
    @test Set(chain.propagator_two_js for chain in all_chains) == Set([SVector(0)])

    @test !isempty(possible_vertex_ls(SpinParity(1, '-'), SpinParity(1, '+'), SpinParity(0, '+')))

    ambiguous_parent = possible_ls_more(SpinParity(1, '+'), SpinParity(1, '+'); jp = SpinParity(2, UndefinedParity))
    definite_plus = possible_ls_more(SpinParity(1, '+'), SpinParity(1, '+'); jp = SpinParity(2, '+'))
    definite_minus = possible_ls_more(SpinParity(1, '+'), SpinParity(1, '+'); jp = SpinParity(2, '-'))
    @test Set(vcat(definite_plus, definite_minus)) ⊆ Set(ambiguous_parent)

    jps = line_spin_parities(topology, weak_system, propagators)
    @test jps[1].p == '+'
    @test jps[root_line_ind(topology)].p == '+'

    undefined_root = CascadeSystem(
        SystemSpinParities(system.quantum, '+', '+', '+'; P0 = UndefinedParity),
        system.masses,
    )
    @test line_spin_parities(topology, undefined_root, propagators)[root_line_ind(topology)].p ==
        UndefinedParity
    undefined_couplings = possible_vertex_couplings(topology, undefined_root, propagators)
    @test !isempty(undefined_couplings[1].second)
    @test !isempty(undefined_couplings[2].second)

    spin_only_resonance = (
        (1, 2) => Propagator(SpinParity(2, '+'), ConstantLineshape(1.0)),
    )
    resonance_couplings = possible_vertex_couplings(topology, weak_system, spin_only_resonance)
    @test isempty(resonance_couplings[1].second)
    @test isempty(resonance_couplings[2].second)
    @test_throws ArgumentError minimal_ls_decay_chain(topology, weak_system, spin_only_resonance)
end
