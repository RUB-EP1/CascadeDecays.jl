struct VertexPrograms{P, L}
    programs::P
    labels::L
end

Base.length(programs::VertexPrograms) = length(programs.programs)
Base.getindex(programs::VertexPrograms, i::Integer) = programs.programs[i]
Base.iterate(programs::VertexPrograms, state...) = iterate(programs.programs, state...)

function Base.show(io::IO, programs::VertexPrograms)
    print(io, "VertexPrograms with ", length(programs), " measurements")
    for (label, program) in zip(programs.labels, programs)
        print(io, "\n  vertex ", label, ": ")
        show(io, program)
    end
    return
end

struct TopologyPrograms{V, A}
    vertex_programs::V
    alignment_paths::A
end

const _REST_FRAME_RTOL = 1.0e-12

"""
    KinematicTask(topologies; reference_topology, wigner_finals, initial_frame)

Reusable kinematic specification: topologies, reference frame for alignment,
which final indices to align, and precompiled IDT programs.
`initial_frame` is stored so evaluation can fall back to `CurrentFrame()` when
the reconstructed parent momentum is numerically already at rest.
"""
struct KinematicTask{Tops, Ref, F, Frame, Progs}
    topologies::Tops
    reference_topology::Ref
    wigner_finals::F
    initial_frame::Frame
    programs::Progs
end

"""
    KinematicPoint

Evaluated kinematic data for one event and one [`KinematicTask`](@ref). It
stores one [`CascadeKinematics`](@ref) object per task topology, plus one
external-axis alignment tuple per topology for the final particles requested by
`task.wigner_finals`.
"""
struct KinematicPoint{Task, Kins, Aligns}
    task::Task
    kinematics::Kins
    alignments::Aligns
end

function KinematicTask(
        topologies::Tuple{Vararg{DecayTopology}};
        reference_topology::DecayTopology = topologies[1],
        wigner_finals::Tuple{Vararg{Int}} = (),
        initial_frame::AbstractInitialFrame = HelicityRootFrame(),
    )
    isempty(topologies) && throw(ArgumentError("KinematicTask requires at least one topology"))
    for final_ind in wigner_finals
        final_ind >= 1 || throw(ArgumentError("wigner_finals entries must be positive final indices"))
    end
    programs = ntuple(length(topologies)) do i
        t = topologies[i]
        vertex_programs = VertexPrograms(
            helicity_angle_programs(t; initial_frame),
            ntuple(v -> _compact_vertex_label(vertex_address(t, v)), nvertices(t)),
        )
        if isempty(wigner_finals)
            return TopologyPrograms(vertex_programs, ())
        end
        alignment_paths = ntuple(length(wigner_finals)) do k
            final_ind = wigner_finals[k]
            final_ind in Base.OneTo(nfinal(t)) ||
                throw(ArgumentError("wigner_finals[$k]=$final_ind outside 1:$(nfinal(t)) for topology $i"))
            (
                helicity_frame_path(reference_topology, final_ind; initial_frame),
                helicity_frame_path(t, final_ind; initial_frame),
            )
        end
        return TopologyPrograms(vertex_programs, alignment_paths)
    end
    return KinematicTask{
        typeof(topologies),
        typeof(reference_topology),
        typeof(wigner_finals),
        typeof(initial_frame),
        typeof(programs),
    }(
        topologies,
        reference_topology,
        wigner_finals,
        initial_frame,
        programs,
    )
end

function _root_momentum(topology::DecayTopology, objs)
    return _sum_objects(objs, _indices_for_line_ind(topology, root_line_ind(topology)))
end

function _line_masses2_from_objects(topology::DecayTopology{Nl}, objs) where {Nl}
    return ntuple(Val(Nl)) do line_ind
        mass(_sum_objects(objs, _indices_for_line_ind(topology, line_ind)))^2
    end
end

function _effectively_at_rest(p; rtol::Real = _REST_FRAME_RTOL)
    scale = max(abs(p.E), one(float(abs(p.E))))
    return p.px^2 + p.py^2 + p.pz^2 <= (rtol * scale)^2
end

function _effective_initial_frame(
        topology::DecayTopology,
        objs,
        initial_frame::AbstractInitialFrame,
    )
    initial_frame isa HelicityRootFrame || return initial_frame
    return _effectively_at_rest(_root_momentum(topology, objs)) ? CurrentFrame() : initial_frame
end

function _topology_slot(task::KinematicTask, topology::DecayTopology)
    idx = findfirst(==(topology), task.topologies)
    idx === nothing &&
        throw(ArgumentError("topology is not part of the kinematic task"))
    return idx
end

"""
    kinematics_at(point, topology)

Return the [`CascadeKinematics`](@ref) stored in `point` for `topology`.
This is a retrieval helper; all kinematic values were computed when the
[`KinematicPoint`](@ref) was built.
"""
kinematics_at(point::KinematicPoint, topology::DecayTopology) =
    point.kinematics[_topology_slot(point.task, topology)]

"""
    alignment_angles_at(point, topology)

Return the relative Wigner alignment angles stored in `point` for `topology`.
Axes follow final-state particle order. Final particles not requested by
`point.task.wigner_finals` carry the identity rotation.
"""
alignment_angles_at(point::KinematicPoint, topology::DecayTopology) =
    point.alignments[_topology_slot(point.task, topology)]

"""
    KinematicPoint(task, objs)

Build one [`KinematicPoint`](@ref) for external four-vectors `objs`. The point
stores one [`CascadeKinematics`](@ref) per task topology plus requested relative
Wigner alignment angles.
"""
function KinematicPoint(task::KinematicTask, objs)
    kinematics = ntuple(length(task.topologies)) do i
        t = task.topologies[i]
        progs = task.programs[i]
        initial_frame = _effective_initial_frame(t, objs, task.initial_frame)
        vertex_programs =
            initial_frame === task.initial_frame ? progs.vertex_programs :
            helicity_angle_programs(t; initial_frame)
        angle_results = map(vertex_programs) do program
            _, result = apply_decay_instruction(program, objs)
            only(values(result))
        end
        CascadeKinematics(_line_masses2_from_objects(t, objs), Tuple(angle_results))
    end
    alignments = ntuple(length(task.topologies)) do i
        t = task.topologies[i]
        progs = task.programs[i]
        Nf = nfinal(t)
        alignments_tuple =
        if isempty(task.wigner_finals)
            ntuple(_ -> _trivial_wigner, Val(Nf))
        else
            initial_frame = _effective_initial_frame(t, objs, task.initial_frame)
            alignment_paths =
                initial_frame === task.initial_frame ? progs.alignment_paths :
                ntuple(length(task.wigner_finals)) do k
                    final_ind = task.wigner_finals[k]
                    (
                        helicity_frame_path(task.reference_topology, final_ind; initial_frame),
                        helicity_frame_path(t, final_ind; initial_frame),
                    )
            end
            ntuple(Val(Nf)) do final_ind
                requested = findfirst(==(final_ind), task.wigner_finals)
                requested === nothing && return _trivial_wigner
                path_ref, path_t = alignment_paths[requested]
                if path_ref == path_t
                    return _trivial_wigner
                end
                cmp = compare_instruction_paths(path_ref, path_t, objs)
                zyz = wigner_zyz(cmp.relative; atol = _WIGNER_DECODE_ATOL)
                return (α = zyz.ϕ, cosβ = cos(zyz.θ), γ = zyz.ψ)
            end
        end
        return SVector{Nf, WignerAngles}(alignments_tuple)
    end
    return KinematicPoint(task, kinematics, alignments)
end

"""
    CascadeKinematics(topology, objs; initial_frame=HelicityRootFrame())

Compute a single-topology [`CascadeKinematics`](@ref) from external four-vectors.
"""
function CascadeKinematics(
        topology::DecayTopology,
        objs;
        initial_frame::AbstractInitialFrame = HelicityRootFrame(),
    )
    task = KinematicTask((topology,); initial_frame = initial_frame)
    return KinematicPoint(task, objs).kinematics[1]
end

"""
    amplitude(chain, system, point::KinematicPoint)

Chain-local helicity amplitude in the topology of `chain`, followed by external
Wigner rotations on the final-state axes listed in `point.task.wigner_finals`
(relative to `point.task.reference_topology`).
"""
function amplitude(
        chain::DecayChain{Nf, Np},
        system::CascadeSystem,
        point::KinematicPoint,
    ) where {Nf, Np}
    x = kinematics_at(point, chain.topology)
    A = _vertex_helicity_amplitude(chain, system, x)
    if !isempty(point.task.wigner_finals)
        A = _apply_external_wigner_rotations(
            A,
            chain,
            system,
            alignment_angles_at(point, chain.topology),
        )
    end
    return A
end

function amplitude(
        chain::DecayChain,
        system::CascadeSystem,
        point::KinematicPoint,
        external_two_λs::SystemSpins,
    )
    A = amplitude(chain, system, point)
    return A[_external_amplitude_indices(chain, system, external_two_λs)...]
end
