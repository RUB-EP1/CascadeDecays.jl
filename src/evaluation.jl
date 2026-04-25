_indices_for_line(topology::DecayTopology, line::Integer) =
    Tuple(final_descendants(topology, line))

_neg_indices(indices::Tuple) = Tuple(-i for i in indices)

function _child_containing_line(topology::DecayTopology, vertex::Integer, line::Integer)
    for child in child_lines(topology, vertex)
        descendants = final_descendants(topology, child)
        if child == line || line in descendants || !isfinalline(topology, line) && !isrootline(topology, line) && line in _subtree_lines(topology, child)
            return child
        end
    end
    throw(ArgumentError("line $line is not below vertex $vertex"))
end

function _subtree_lines(topology::DecayTopology, line::Integer)
    lines = Int[Int(line)]
    vertex = consumed_by(topology, line)
    vertex === nothing && return lines
    for child in child_lines(topology, vertex)
        append!(lines, _subtree_lines(topology, child))
    end
    return lines
end

function _root_vertex(topology::DecayTopology)
    vertex = consumed_by(topology, rootline(topology))
    vertex === nothing && throw(ArgumentError("topology root line is not consumed by any vertex"))
    return vertex
end

abstract type AbstractInitialFrame end

"""
    HelicityRootFrame()

Start generated helicity-angle programs by transforming to the root/system
helicity frame. This is the default for fully general four-vectors.
"""
struct HelicityRootFrame <: AbstractInitialFrame end

"""
    CurrentFrame()

Start generated helicity-angle programs in the current axes. Use this when the
input four-vectors are already expressed in the intended system frame.
"""
struct CurrentFrame <: AbstractInitialFrame end

_initial_frame_program(topology::DecayTopology, ::HelicityRootFrame) =
    (ToHelicityFrame(_indices_for_line(topology, rootline(topology))),)
_initial_frame_program(::DecayTopology, ::CurrentFrame) = ()

"""
    helicity_angle_program(topology, vertex; initial_frame=HelicityRootFrame())

Build an `InstructionalDecayTrees.jl` instruction program that measures the
local `(cosθ, ϕ)` angle for `vertex` in helicity convention.
"""
function helicity_angle_program(
    topology::DecayTopology,
    target_vertex::Integer;
    initial_frame::AbstractInitialFrame = HelicityRootFrame(),
)
    return helicity_angle_program(topology, Val(Int(target_vertex)); initial_frame)
end

function helicity_angle_program(
    topology::DecayTopology,
    ::Val{target_vertex};
    initial_frame::AbstractInitialFrame = HelicityRootFrame(),
) where {target_vertex}
    target_vertex isa Integer ||
        throw(ArgumentError("vertex $target_vertex is outside 1:$(nvertices(topology))"))
    target_vertex in Base.OneTo(nvertices(topology)) ||
        throw(ArgumentError("vertex $target_vertex is outside 1:$(nvertices(topology))"))

    target_parent = incoming_line(topology, target_vertex)
    program = _initial_frame_program(topology, initial_frame)

    current_vertex = _root_vertex(topology)

    while true
        parent = incoming_line(topology, current_vertex)
        children = child_lines(topology, current_vertex)
        if parent == target_parent
            return (
                program...,
                MeasureCosThetaPhi(Symbol(:v, target_vertex), _indices_for_line(topology, children[1])),
            )
        end
        next_child = _child_containing_line(topology, current_vertex, target_parent)
        program = (program..., ToHelicityFrame(_indices_for_line(topology, next_child)))
        current_vertex = consumed_by(topology, next_child)
        current_vertex === nothing &&
            throw(ArgumentError("target vertex $target_vertex is not reachable from the root"))
    end
end

"""
    helicity_angle_programs(topology; initial_frame=HelicityRootFrame())

Return one helicity-angle measurement program per topology vertex.
"""
helicity_angle_programs(
    topology::DecayTopology;
    initial_frame::AbstractInitialFrame = HelicityRootFrame(),
) =
    ntuple(v -> helicity_angle_program(topology, Val(v); initial_frame), nvertices(topology))

function _sum_objects(objs, indices::Tuple)
    return reduce(+, (objs[i] for i in indices))
end

"""
    cascade_kinematics(topology, system, objs; initial_frame=HelicityRootFrame())

Compute a `CascadeKinematics` input from final-state four-vectors `objs`.
Internal invariant masses are derived from final-state descendants, and vertex
angles are computed with generated helicity-angle programs.
"""
function cascade_kinematics(
    topology::DecayTopology,
    system::CascadeSystem,
    objs;
    initial_frame::AbstractInitialFrame = HelicityRootFrame(),
)
    internal_masses2 = Tuple(mass(_sum_objects(objs, _indices_for_line(topology, line)))^2 for line in internal_lines(topology))
    programs = helicity_angle_programs(topology; initial_frame)
    angle_results = map(programs) do program
        _, result = apply_decay_instruction(program, objs)
        only(values(result))
    end
    return CascadeKinematics(topology, system; internal_masses2, vertex_angles = Tuple(angle_results))
end

function routed_vertex_amplitude(vertex, masses2, helicities, spins, angles)
    throw(MethodError(routed_vertex_amplitude, (vertex, masses2, helicities, spins, angles)))
end

function routed_vertex_amplitude(
    vertex::ThreeBodyDecays.VertexFunction,
    masses2,
    helicities,
    spins,
    angles,
)
    two_j0, _, _ = spins
    two_λ0, two_λ1, two_λ2 = helicities
    recoupling =
        _particle_two_phase(spins[3], two_λ2) *
        ThreeBodyDecays.amplitude(vertex.h, (two_λ1, two_λ2), spins)
    formfactor = vertex.ff(masses2...)
    two_Δλ = two_λ1 - two_λ2
    rotation = conj(ThreeBodyDecays.wignerD_doublearg(two_j0, two_λ0, two_Δλ, angles.ϕ, angles.cosθ, 0))
    return rotation * recoupling * formfactor
end

function _particle_two_phase(two_j2::Integer, two_λ2::Integer)
    exponent_num = two_j2 - two_λ2
    iseven(exponent_num) ||
        throw(ArgumentError("particle-2 phase requires two_j2 - two_λ2 to be even"))
    return isodd(div(exponent_num, 2)) ? -1 : 1
end

function routed_vertex_amplitude(chain::DecayChain, system::CascadeSystem, x::CascadeKinematics, two_λs, vertex::Integer)
    masses2 = vertex_masses2(chain, x, vertex)
    helicities = vertex_helicities(chain, two_λs, vertex)
    spins = vertex_spins(chain, system, vertex)
    angles = vertex_angles(x, vertex)
    return routed_vertex_amplitude(chain.vertices[vertex], masses2, helicities, spins, angles)
end

function routed_propagator_product(chain::DecayChain, x::CascadeKinematics)
    return prod(zip(chain.propagators, propagating_lines(chain))) do (propagator, line)
        propagator(line_invariant(x, line))
    end
end

function _helicity_axis(two_j::Integer)
    return (-Int(two_j)):2:Int(two_j)
end

function _full_helicity_assignment(chain::DecayChain, external_two_λs, internal_two_λs)
    external_tuple = Tuple(Int(two_λ) for two_λ in external_two_λs)
    length(external_tuple) == nfinal(chain) + 1 ||
        throw(ArgumentError("external_two_λs must be given as `(finals..., root)`"))
    internal_tuple = Tuple(Int(two_λ) for two_λ in internal_two_λs)
    length(internal_tuple) == length(propagating_lines(chain)) ||
        throw(ArgumentError("internal helicity assignment does not match the number of internal lines"))
    values = MVector{nlines(chain),Int}(undef)
    for (i, line) in pairs(finallines(chain))
        values[line] = external_tuple[i]
    end
    for (i, line) in pairs(propagating_lines(chain))
        values[line] = internal_tuple[i]
    end
    values[rootline(chain)] = external_tuple[end]
    return SVector(values)
end

function _internal_helicity_assignments(chain::DecayChain)
    axes = Tuple(_helicity_axis(two_j) for two_j in chain.propagator_two_js)
    isempty(axes) && return ((),)
    return Iterators.product(axes...)
end

"""
    amplitude(chain, system, x, external_two_λs)

Route masses, angles, and spins through the cascade and sum over all internal
helicity assignments. `external_two_λs` is supplied as `(finals..., root)`.
"""
function amplitude(chain::DecayChain, system::CascadeSystem, x::CascadeKinematics, external_two_λs)
    P = routed_propagator_product(chain, x)
    return sum(_internal_helicity_assignments(chain)) do internal_two_λs
        two_λs = _full_helicity_assignment(chain, external_two_λs, internal_two_λs)
        V = prod(v -> routed_vertex_amplitude(chain, system, x, two_λs, v), 1:nvertices(chain))
        V * P
    end
end
