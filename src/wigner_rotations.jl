const _WIGNER_DECODE_ATOL = 1e-10

"""
    WignerAngles

Standard ZYZ Wigner angles `(α, cosβ, γ)` for an external helicity axis.
Identity rotation: `(α=0, cosβ=1, γ=0)`.
"""
const WignerAngles = NamedTuple{(:α, :cosβ, :γ), NTuple{3, Float64}}
const _trivial_wigner = (α = 0.0, cosβ = 1.0, γ = 0.0)

function _helicity_range(two_j::Integer)
    return (-two_j):2:two_j
end

function wigner_d_zyz_matrix(two_j::Integer, α::Real, cosβ::Real, γ::Real)
    rows = _helicity_range(two_j)
    return wignerD_doublearg.(two_j, rows, transpose(rows), α, cosβ, γ)
end

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
    helicity_frame_paths(topology; initial_frame=HelicityRootFrame())

Return one helicity-frame path per final-state particle in canonical order.
"""
helicity_frame_paths(
    topology::DecayTopology;
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
) = ntuple(i -> helicity_frame_path(topology, i; initial_frame), nfinal(topology))

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

function _root_relative_wigner_angles(
    reference_topology::DecayTopology,
    topology::DecayTopology,
    objs;
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
    T::Type{<:Real}=Float64,
)
    path_ref = _initial_frame_program(reference_topology, initial_frame)
    path_other = _initial_frame_program(topology, initial_frame)
    path_ref == path_other && return _trivial_alignment_rotation
    cmp = compare_instruction_paths(path_ref, path_other, objs; T = T)
    return wigner_zyz(cmp.relative; atol = _WIGNER_DECODE_ATOL)
end

"""
    external_wigner_angles(reference_topology, topology, objs; initial_frame=HelicityRootFrame())

Precompute relative Wigner angles for every external axis of `topology` in
[`external_line_inds`](@ref) order: final-state particles followed by the root.
"""
function external_wigner_angles(
    reference_topology::DecayTopology,
    topology::DecayTopology,
    objs;
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
)
    final_angles = ntuple(nfinal(topology)) do i
        relative_wigner_angles(
            reference_topology,
            topology,
            i,
            objs;
            initial_frame,
        )
    end
    root_angles = _root_relative_wigner_angles(
        reference_topology,
        topology,
        objs;
        initial_frame,
    )
    return (final_angles..., root_angles)
end

_wigner_d_matrix(two_j::Integer, angles::WignerAngles) =
    conj.(wigner_d_zyz_matrix(two_j, angles.α, angles.cosβ, angles.γ))

function _external_wigner_matrices(
    chain::DecayChain,
    system::CascadeSystem,
    angles::SVector{Nf,WignerAngles},
) where {Nf}
    final_lines = final_line_inds(chain.topology)
    length(angles) == length(final_lines) ||
        throw(ArgumentError("angles must have one entry per final-state particle"))
    ext_lines = external_line_inds(chain.topology)
    two_js = line_two_js(chain, system)
    return ntuple(length(ext_lines)) do i
        line = ext_lines[i]
        wigner_angles =
            line == root_line_ind(chain.topology) ? _trivial_wigner :
            angles[findfirst(==(line), final_lines)]
        two_j = two_js[line]
        two_j == 0 ? ones(Float64, 1, 1) : _wigner_d_matrix(two_j, wigner_angles)
    end
end

# Generated `@tullio` contractions for external Wigner rotations.
# External axis order is (finals..., root); each axis uses ``F[m',\\ldots] \\cdot D[m',m]``.
for Ne in 1:8
    ext_syms = Tuple(Symbol(:_λ, k) for k in 1:Ne)
    extp_syms = Tuple(Symbol(:_λp, k) for k in 1:Ne)
    D_refs = [:($(Symbol(:D, k))[$(extp_syms[k]), $(ext_syms[k])]) for k in 1:Ne]
    D_unpack = map(k -> :($(Symbol(:D, k)) = Ds[$k]), 1:Ne)
    F_ref = Expr(:ref, :F, extp_syms...)
    A_ref = Expr(:ref, :A, ext_syms...)
    tullio_rhs = reduce((a, b) -> :($a * $b), D_refs; init = F_ref)
    @eval function _apply_external_wigner_rotations(
        F::AbstractArray{Ta,$Ne},
        Ds::Tuple,
    ) where {Ta}
        length(Ds) == $Ne || throw(ArgumentError("expected $Ne Wigner matrices"))
        $(D_unpack...)
        A = zeros(Ta, size(F)...)
        @tullio $A_ref := $tullio_rhs
        return A
    end
end

function _apply_external_wigner_rotations(F::AbstractArray{Ta,Ne}, Ds::Tuple) where {Ta,Ne}
    throw(ArgumentError("unsupported external axis count Ne=$Ne for Wigner contraction"))
end

function apply_external_wigner_rotations(
    F::AbstractArray,
    chain::DecayChain,
    system::CascadeSystem,
    angles::SVector{Nf,WignerAngles},
) where {Nf}
    Ds = _external_wigner_matrices(chain, system, angles)
    all(D -> size(D) == (1, 1), Ds) && return F
    return _apply_external_wigner_rotations(F, Ds)
end
