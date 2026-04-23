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

"""
    helicity_angle_program(topology, vertex)

Build an `InstructionalDecayTrees.jl` instruction program that measures the
local `(cosθ, ϕ)` angle for `vertex` in helicity convention.
"""
function helicity_angle_program(topology::DecayTopology, target_vertex::Integer)
    target_vertex in Base.OneTo(nvertices(topology)) ||
        throw(ArgumentError("vertex $target_vertex is outside 1:$(nvertices(topology))"))

    target_parent = incoming_line(topology, target_vertex)
    program = Any[ToHelicityFrame(_indices_for_line(topology, rootline(topology)))]

    current_vertex = _root_vertex(topology)
    root_children = child_lines(topology, current_vertex)
    push!(
        program,
        PlaneAlign(
            _neg_indices(_indices_for_line(topology, root_children[2])),
            _indices_for_line(topology, root_children[1]),
        ),
    )

    while true
        parent = incoming_line(topology, current_vertex)
        children = child_lines(topology, current_vertex)
        if parent == target_parent
            push!(program, MeasureCosThetaPhi(Symbol(:v, target_vertex), _indices_for_line(topology, children[1])))
            return Tuple(program)
        end
        next_child = _child_containing_line(topology, current_vertex, target_parent)
        push!(program, ToHelicityFrame(_indices_for_line(topology, next_child)))
        current_vertex = consumed_by(topology, next_child)
        current_vertex === nothing &&
            throw(ArgumentError("target vertex $target_vertex is not reachable from the root"))
    end
end

"""
    helicity_angle_programs(topology)

Return one helicity-angle measurement program per topology vertex.
"""
helicity_angle_programs(topology::DecayTopology) =
    ntuple(v -> helicity_angle_program(topology, v), nvertices(topology))

function _sum_objects(objs, indices::Tuple)
    return reduce(+, (objs[i] for i in indices))
end

"""
    cascade_kinematics(topology, system, objs)

Compute a `CascadeKinematics` input from final-state four-vectors `objs`.
Internal invariant masses are derived from final-state descendants, and vertex
angles are computed with generated helicity-angle programs.
"""
function cascade_kinematics(topology::DecayTopology, system::CascadeSystem, objs)
    internal_masses2 = Tuple(mass(_sum_objects(objs, _indices_for_line(topology, line)))^2 for line in internal_lines(topology))
    programs = helicity_angle_programs(topology)
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
    recoupling = ThreeBodyDecays.amplitude(vertex.h, (two_λ1, two_λ2), spins)
    formfactor = vertex.ff(masses2...)
    two_Δλ = two_λ1 - two_λ2
    rotation = conj(ThreeBodyDecays.wignerD_doublearg(two_j0, two_λ0, two_Δλ, angles.ϕ, angles.cosθ, 0))
    return rotation * recoupling * formfactor
end

function routed_vertex_amplitude(chain::DecayChain, system::CascadeSystem, x::CascadeKinematics, two_λs, vertex::Integer)
    masses2 = vertex_masses2(chain, x, vertex)
    helicities = vertex_helicities(chain, two_λs, vertex)
    spins = vertex_spins(chain.topology, system, vertex)
    angles = vertex_angles(x, vertex)
    return routed_vertex_amplitude(chain.vertices[vertex], masses2, helicities, spins, angles)
end

function routed_propagator_product(chain::DecayChain, x::CascadeKinematics)
    return prod(zip(chain.propagators, propagating_lines(chain))) do (propagator, line)
        propagator(line_invariant(x, line))
    end
end

"""
    amplitude(chain, system, x, two_λs)

Route graph-indexed masses, angles, spins, and helicities to static vertex and
propagator payloads and return the product for one complete helicity assignment.
"""
function amplitude(chain::DecayChain, system::CascadeSystem, x::CascadeKinematics, two_λs)
    V = prod(v -> routed_vertex_amplitude(chain, system, x, two_λs, v), 1:nvertices(chain))
    P = routed_propagator_product(chain, x)
    return V * P
end
