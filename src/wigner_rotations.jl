const _WIGNER_DECODE_ATOL = 1e-10

"""
    WignerAngles

Standard ZYZ Wigner angles `(α, cosβ, γ)` for an external helicity axis.
Identity rotation: `(α=0, cosβ=1, γ=0)`.
"""
const WignerAngles = NamedTuple{(:α, :cosβ, :γ), NTuple{3, Float64}}
const _trivial_wigner = (α = 0.0, cosβ = 1.0, γ = 0.0)

function _path_steps_to_line(topology::DecayTopology, line_ind::Integer)
    line_ind == root_line_ind(topology) && return ()
    path = _find_path_to_line(topology, root_line_ind(topology), Int(line_ind))
    path === nothing &&
        throw(ArgumentError("line_ind $line_ind is not reachable from the root"))
    return path
end

function _find_path_to_line(topology::DecayTopology, current_line::Int, target::Int)
    current_line == target && return ()
    vertex_ind = consumed_by(topology, current_line)
    vertex_ind === nothing && return nothing
    for child in child_line_inds(topology, vertex_ind)
        path = _find_path_to_line(topology, child, target)
        path === nothing && continue
        return ((vertex_ind = vertex_ind, child_line = child), path...)
    end
    return nothing
end

function _helicity_step_instruction(topology::DecayTopology, vertex_ind::Integer, child_line::Integer)
    indices = _indices_for_line_ind(topology, child_line)
    _, child1, child2 = vertex_line_inds(topology, vertex_ind)
    return child_line == child2 ? ToHelicityFrameParticle2(indices) : ToHelicityFrame(indices)
end

"""
    helicity_frame_path(topology, particle_index; initial_frame=HelicityRootFrame())

Build an `InstructionalDecayTrees.jl` instruction path that defines the helicity
quantization frame for final-state particle `particle_index`. A final particle
is represented internally by a topology line; this method keeps that line id out
of the user-facing API.
"""
function helicity_frame_path(
    topology::DecayTopology,
    particle_index::Integer;
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
)
    particle_index in Base.OneTo(nfinal(topology)) ||
        throw(ArgumentError("particle_index $particle_index is outside 1:$(nfinal(topology))"))
    line_ind = final_line_inds(topology)[particle_index]
    _require_line_ind(topology, line_ind)
    program = _initial_frame_program(topology, initial_frame)
    for step in _path_steps_to_line(topology, line_ind)
        program = (
            program...,
            _helicity_step_instruction(topology, step.vertex_ind, step.child_line),
        )
    end
    return program
end

"""
    relative_wigner_angles(reference_topology, topology, particle_index, objs; T=Float64)

Compare helicity-frame instruction paths for final-state particle
`particle_index` and return ZYZ Wigner angles `(ϕ, θ, ψ)` of the rotation from
the reference topology frame to the `topology` frame.
"""
function relative_wigner_angles(
    reference_topology::DecayTopology,
    topology::DecayTopology,
    particle_index::Integer,
    objs;
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
    T::Type{<:Real}=Float64,
)
    path_ref = helicity_frame_path(reference_topology, particle_index; initial_frame)
    path_other = helicity_frame_path(topology, particle_index; initial_frame)
    path_ref == path_other && return _trivial_alignment_rotation
    cmp = compare_instruction_paths(path_ref, path_other, objs; T = T)
    return wigner_zyz(cmp.relative; atol = _WIGNER_DECODE_ATOL)
end
