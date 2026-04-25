using ThreeBodyDecays: ThreeBodyMasses, aligned_four_vectors, cosθij, x2σs

function _fourvector_from_tuple(p)
    return FourVector(Float64(p[1]), Float64(p[2]), Float64(p[3]); E = Float64(p[4]))
end

function _threebody_reference_point()
    ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0 = 7.7)
    σs = x2σs([0.5, 0.3], ms; k = 3)
    topology = DecayTopology(((1, 2), 3))
    system = CascadeSystem((0, 0, 0, 0), (ms[1]^2, ms[2]^2, ms[3]^2, ms.m0^2))
    return (; ms, σs, topology, system)
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
    generated_x = cascade_kinematics(topology, system, objs)

    @test line_invariant(topology, generated_x, (1, 2)) ≈ σs[3]
    @test vertex_angles(topology, generated_x, ((1, 2), 3)).cosθ ≈ 1.0
    @test vertex_angles(topology, generated_x, ((1, 2), 3)).ϕ ≈ 0.0
    @test vertex_angles(topology, generated_x, (1, 2)).cosθ ≈ expected_inner_cosθ
    @test vertex_angles(topology, generated_x, (1, 2)).ϕ ≈ 0.0
end

@testset "Rotated aligned root kinematics diagnostic" begin
    (; ms, σs, topology, system) = _threebody_reference_point()
    objs = _fourvector_from_tuple.(aligned_four_vectors(σs, ms; k = 3))

    θ = deg2rad(10.0)
    ϕ = deg2rad(20.0)
    rotated_objs = map(p -> p |> Ry(θ) |> Rz(ϕ), objs)
    generated_x = cascade_kinematics(topology, system, rotated_objs; initial_frame = CurrentFrame())
    root_angles = vertex_angles(topology, generated_x, ((1, 2), 3))

    @test root_angles.cosθ ≈ cos(θ)
    @test root_angles.ϕ ≈ ϕ
end
