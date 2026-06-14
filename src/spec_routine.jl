struct TopologyPrograms
    vertex_programs::Tuple
    alignment_paths::Tuple
end

const _REST_FRAME_RTOL = 1e-12

"""
    KinematicTask(topologies; reference_topology, wigner_finals, initial_frame)

Reusable kinematic specification: topologies, reference frame for alignment,
which final indices to measure at evaluate, and precompiled IDT programs.
`initial_frame` is stored so evaluation can fall back to `CurrentFrame()` when
the reconstructed parent momentum is numerically already at rest.
"""
struct KinematicTask{Tops,Ref,F,Frame,Progs}
    topologies::Tops
    reference_topology::Ref
    wigner_finals::F
    initial_frame::Frame
    programs::Progs
end

struct KinematicPoint{Task,Kins,Aligns}
    task::Task
    kinematics::Kins
    alignments::Aligns
end

function KinematicTask(
    topologies::Tuple{Vararg{DecayTopology}};
    reference_topology::DecayTopology=topologies[1],
    wigner_finals::Tuple{Vararg{Int}}=(),
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
)
    isempty(topologies) && throw(ArgumentError("KinematicTask requires at least one topology"))
    for final_ind in wigner_finals
        final_ind >= 1 || throw(ArgumentError("wigner_finals entries must be positive final indices"))
    end
    programs = ntuple(length(topologies)) do i
        t = topologies[i]
        vertex_programs = helicity_angle_programs(t; initial_frame)
        if isempty(wigner_finals)
            return TopologyPrograms(vertex_programs, ())
        end
        alignment_paths = ntuple(length(wigner_finals)) do k
            final_ind = wigner_finals[k]
            final_ind in Base.OneTo(nfinal(t)) ||
                throw(ArgumentError("wigner_finals[$k]=$final_ind outside 1:$(nfinal(t)) for topology $i"))
            line = final_line_inds(t)[final_ind]
            (
                helicity_frame_path(reference_topology, line; initial_frame),
                helicity_frame_path(t, line; initial_frame),
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

function _effectively_at_rest(p; rtol::Real=_REST_FRAME_RTOL)
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

kinematics_at(point::KinematicPoint, topology::DecayTopology) =
    point.kinematics[_topology_slot(point.task, topology)]

alignment_angles_at(point::KinematicPoint, topology::DecayTopology) =
    point.alignments[_topology_slot(point.task, topology)]

"""
    evaluate(task, objs, system)

Evaluate one kinematic point for external four-vectors `objs` and return per-topology
kinematics plus requested relative Wigner alignment angles.
"""
function evaluate(task::KinematicTask, objs, system::CascadeSystem)
    kinematics = ntuple(length(task.topologies)) do i
        t = task.topologies[i]
        progs = task.programs[i]
        initial_frame = _effective_initial_frame(t, objs, task.initial_frame)
        vertex_programs =
            initial_frame === task.initial_frame ? progs.vertex_programs :
            helicity_angle_programs(t; initial_frame)
        internal_masses2 = Tuple(
            mass(_sum_objects(objs, _indices_for_line_ind(t, line_ind)))^2 for
            line_ind in internal_line_inds(t)
        )
        angle_results = map(vertex_programs) do program
            _, result = apply_decay_instruction(program, objs)
            only(values(result))
        end
        CascadeKinematics(
            t,
            system;
            internal_masses2,
            vertex_angles = Tuple(angle_results),
        )
    end
    alignments = ntuple(length(task.topologies)) do i
        t = task.topologies[i]
        progs = task.programs[i]
        ext_lines = external_line_inds(t)
        Ne = length(ext_lines)
        alignments = fill(_trivial_wigner, Ne)
        if !isempty(task.wigner_finals)
            initial_frame = _effective_initial_frame(t, objs, task.initial_frame)
            alignment_paths =
                initial_frame === task.initial_frame ? progs.alignment_paths :
                ntuple(length(task.wigner_finals)) do k
                    final_ind = task.wigner_finals[k]
                    line = final_line_inds(t)[final_ind]
                    (
                        helicity_frame_path(task.reference_topology, line; initial_frame),
                        helicity_frame_path(t, line; initial_frame),
                    )
                end
            for (k, final_ind) in enumerate(task.wigner_finals)
                path_ref, path_t = alignment_paths[k]
                angles =
                    path_ref == path_t ? _trivial_wigner :
                    begin
                        cmp = compare_instruction_paths(path_ref, path_t, objs)
                        zyz = wigner_zyz(cmp.relative; atol = _WIGNER_DECODE_ATOL)
                        (α = zyz.ϕ, cosβ = cos(zyz.θ), γ = zyz.ψ)
                    end
                line = final_line_inds(t)[final_ind]
                axis = findfirst(==(line), ext_lines)
                axis === nothing &&
                    throw(ArgumentError("final line $line is not an external line of topology"))
                alignments[axis] = angles
            end
        end
        return SVector{Ne,WignerAngles}(alignments)
    end
    return KinematicPoint(task, kinematics, alignments)
end

"""
    cascade_kinematics(topology, system, objs; initial_frame=HelicityRootFrame())

Compute a single-topology [`CascadeKinematics`](@ref) from external four-vectors.
"""
function cascade_kinematics(
    topology::DecayTopology,
    system::CascadeSystem,
    objs;
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
)
    task = KinematicTask((topology,); initial_frame = initial_frame)
    return evaluate(task, objs, system).kinematics[1]
end
