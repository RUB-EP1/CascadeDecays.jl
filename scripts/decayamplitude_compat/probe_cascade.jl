using CascadeDecays
using FourVectors
using TOML
using ThreeBodyDecays: RecouplingLS

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FIXTURES = joinpath(ROOT, "docs", "fixtures", "decayamplitude_compat")

event = TOML.parsefile(joinpath(FIXTURES, "event.toml"))
reference = TOML.parsefile(joinpath(FIXTURES, "reference.toml"))

momenta = Tuple(
    FourVector(p["px"], p["py"], p["pz"]; E = p["E"])
    for p in sort(event["momenta"]; by = p -> p["particle"])
)

topology = DecayTopology(((1, 2), 3))
kinematics = DecayChainKinematics(topology, momenta; initial_frame = CurrentFrame())
spins = SystemSpins(0, 1, 1; two_h0 = 2)
chain = DecayChain(
    topology,
    spins;
    propagators = (((1, 2) => Propagator(1, ConstantLineshape(1.0))),),
    vertices = (
        ((1, 2), 3) => Vertex(RecouplingLS((0, 2))),
        (1, 2) => Vertex(RecouplingLS((0, 1))),
    ),
)

A = amplitude(chain, kinematics)

alternate_topology = DecayTopology(((1, 3), 2))
alternate_chain = DecayChain(
    alternate_topology,
    spins;
    propagators = (((1, 3) => Propagator(1, ConstantLineshape(1.0))),),
    vertices = (
        ((1, 3), 2) => Vertex(RecouplingLS((0, 2))),
        (1, 3) => Vertex(RecouplingLS((0, 1))),
    ),
)
task = KinematicTask(
    (topology, alternate_topology);
    reference_topology = topology,
    wigner_finals = (2, 3),
    initial_frame = CurrentFrame(),
)
point = KinematicPoint(task, momenta)
aligned_A = amplitude(alternate_chain, point)

function amplitude_entry(A, row)
    i0 = div(row["two_h0"] + 2, 2) + 1
    i1 = 1
    i2 = div(row["two_h2"] + 1, 2) + 1
    i3 = div(row["two_h3"] + 1, 2) + 1
    return A[i1, i2, i3, i0]
end

angle_root = vertex_angles(kinematics, topology, ((1, 2), 3))
angle_isobar = vertex_angles(kinematics, topology, (1, 2))
kinematic_errors = vcat(
    abs.(collect(kinematics.line_masses2) .- reference["line_masses2"]),
    [
        abs(angle_root.ϕ - reference["root_phi"]),
        abs(angle_root.cosθ - reference["root_cos_theta"]),
        abs(angle_isobar.ϕ - reference["isobar_phi"]),
        abs(angle_isobar.cosθ - reference["isobar_cos_theta"]),
    ],
)
amplitude_errors = [
    abs(amplitude_entry(A, row) - complex(row["re"], row["im"]))
    for row in reference["amplitudes"]
]
aligned_amplitude_errors = [
    abs(amplitude_entry(aligned_A, row) - complex(row["re"], row["im"]))
    for row in reference["aligned_amplitudes"]
]
alignment_values = alignment_angles_at(point, alternate_topology)
alignment_errors = Float64[]
for row in reference["alignments"]
    value = alignment_values[row["particle"]]
    append!(alignment_errors, [
        abs(value.α - row["alpha"]),
        abs(value.cosβ - row["cos_beta"]),
        abs(value.γ - row["gamma"]),
    ])
end

representative_terms = Dict{Int, NamedTuple}()
for two_λR in (-1, 1)
    two_λs = [0, -1, 1, two_λR, 0]
    root_coupling = CascadeDecays._vertex_coupling_value(
        chain.vertices[1],
        vertex_masses2(chain, kinematics, 1),
        two_λR,
        1,
        CascadeDecays.vertex_spins(chain, 1),
    )
    isobar_coupling = CascadeDecays._vertex_coupling_value(
        chain.vertices[2],
        vertex_masses2(chain, kinematics, 2),
        0,
        -1,
        CascadeDecays.vertex_spins(chain, 2),
    )
    root_amplitude = routed_vertex_amplitude(chain, kinematics, two_λs, 1)
    isobar_amplitude = routed_vertex_amplitude(chain, kinematics, two_λs, 2)
    representative_terms[two_λR] = (;
        root_coupling,
        isobar_coupling,
        root_amplitude,
        isobar_amplitude,
        product = root_amplitude * isobar_amplitude,
    )
end

term_errors = Float64[]
for row in reference["terms"]
    term = representative_terms[row["two_lambda_internal"]]
    for name in (:root_coupling, :isobar_coupling, :root_amplitude, :isobar_amplitude, :product)
        re = row[string(name, "_re")]
        im = row[string(name, "_im")]
        push!(term_errors, abs(getproperty(term, name) - complex(re, im)))
    end
end

max_kinematic_error = maximum(kinematic_errors)
max_term_error = maximum(term_errors)
max_amplitude_error = maximum(amplitude_errors)
max_alignment_error = maximum(alignment_errors)
max_aligned_amplitude_error = maximum(aligned_amplitude_errors)

@assert max_kinematic_error < 1.0e-14
@assert max_term_error < 1.0e-14
@assert max_amplitude_error < 1.0e-14
@assert max_alignment_error < 1.0e-14
@assert max_aligned_amplitude_error < 1.0e-14

println("CascadeDecays one-event probe")
println("  topology: ", bracket_notation(topology))
println("  max kinematic error: ", max_kinematic_error)
println("  max computation-graph term error: ", max_term_error)
println("  max external-amplitude error: ", max_amplitude_error)
println("  max alignment-angle error: ", max_alignment_error)
println("  max aligned-amplitude error: ", max_aligned_amplitude_error)
