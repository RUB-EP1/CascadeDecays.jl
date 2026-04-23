"""
    CascadeKinematics(line_masses2, vertex_angles)

Runtime kinematic input point. `line_masses2` is indexed by line id and
`vertex_angles` is indexed by vertex id.
"""
struct CascadeKinematics{Nl,Nv,T,A}
    line_masses2::SVector{Nl,T}
    vertex_angles::SVector{Nv,A}
end

_has_vertex_angle_schema(angle) = hasproperty(angle, :cosθ) && hasproperty(angle, :ϕ)

function CascadeKinematics(line_masses2, vertex_angles)
    mass_tuple = Tuple(line_masses2)
    angle_tuple = Tuple(vertex_angles)
    all(_has_vertex_angle_schema, angle_tuple) ||
        throw(ArgumentError("each vertex angle must provide `cosθ` and `ϕ` fields"))
    T = promote_type(map(typeof, mass_tuple)...)
    A = promote_type(map(typeof, angle_tuple)...)
    return CascadeKinematics{length(mass_tuple),length(angle_tuple),T,A}(
        SVector{length(mass_tuple),T}(mass_tuple),
        SVector{length(angle_tuple),A}(angle_tuple),
    )
end

function CascadeKinematics(
    topology::DecayTopology,
    system::CascadeSystem;
    internal_masses2,
    vertex_angles,
)
    masses = line_masses2(topology, system, internal_masses2)
    length(vertex_angles) == nvertices(topology) ||
        throw(ArgumentError("vertex_angles must have one entry per topology vertex"))
    return CascadeKinematics(masses, vertex_angles)
end

function _require_kinematic_line(x::CascadeKinematics, line::Integer)
    line in Base.OneTo(length(x.line_masses2)) ||
        throw(ArgumentError("line $line is outside 1:$(length(x.line_masses2))"))
    return nothing
end

function _require_kinematic_vertex(x::CascadeKinematics, vertex::Integer)
    vertex in Base.OneTo(length(x.vertex_angles)) ||
        throw(ArgumentError("vertex $vertex is outside 1:$(length(x.vertex_angles))"))
    return nothing
end

function _check_kinematics(topology::DecayTopology, x::CascadeKinematics{Nl,Nv}) where {Nl,Nv}
    Nl == nlines(topology) ||
        throw(ArgumentError("kinematics line_masses2 must have one entry per topology line"))
    Nv == nvertices(topology) ||
        throw(ArgumentError("kinematics vertex_angles must have one entry per topology vertex"))
    return nothing
end

function line_invariant(x::CascadeKinematics, line::Integer)
    _require_kinematic_line(x, line)
    return x.line_masses2[line]
end

function vertex_masses2(topology::DecayTopology, x::CascadeKinematics, vertex::Integer)
    _check_kinematics(topology, x)
    l0, l1, l2 = vertex_lines(topology, vertex)
    return (x.line_masses2[l0], x.line_masses2[l1], x.line_masses2[l2])
end

function vertex_helicities(topology::DecayTopology, two_λs, vertex::Integer)
    l0, l1, l2 = vertex_lines(topology, vertex)
    length(two_λs) == nlines(topology) ||
        throw(ArgumentError("two_λs must have one entry per topology line"))
    return (two_λs[l0], two_λs[l1], two_λs[l2])
end

function vertex_spins(topology::DecayTopology, system::CascadeSystem, vertex::Integer)
    _check_system(topology, system)
    l0, l1, l2 = vertex_lines(topology, vertex)
    return (system.two_js[l0], system.two_js[l1], system.two_js[l2])
end

function vertex_angles(x::CascadeKinematics, vertex::Integer)
    _require_kinematic_vertex(x, vertex)
    return x.vertex_angles[vertex]
end
