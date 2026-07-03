"""
    DecayChainKinematics{Nl, Nv, T, A}

Runtime kinematic input point for one decay topology and one event.

The object stores line invariant masses squared and local helicity angles in
the internal indexing of a [`DecayTopology`](@ref). Prefer the retrieval helpers
[`line_invariant`](@ref), [`vertex_masses2`](@ref), and [`vertex_angles`](@ref),
which accept line or vertex ids and, with a topology, bracket addresses.

# Fields

- `line_masses2::SVector{Nl, T}`: invariant masses squared indexed by line id.
- `vertex_angles::SVector{Nv, A}`: helicity angles indexed by vertex id.
"""
struct DecayChainKinematics{Nl, Nv, T, A}
    line_masses2::SVector{Nl, T}
    vertex_angles::SVector{Nv, A}
end

_has_vertex_angle_schema(angle) = hasproperty(angle, :cosθ) && hasproperty(angle, :ϕ)

function _decay_chain_kinematics_from_values(line_masses2, vertex_angles)
    mass_tuple = Tuple(line_masses2)
    angle_tuple = Tuple(vertex_angles)
    all(_has_vertex_angle_schema, angle_tuple) ||
        throw(ArgumentError("each vertex angle must provide `cosθ` and `ϕ` fields"))
    T = promote_type(map(typeof, mass_tuple)...)
    A = promote_type(map(typeof, angle_tuple)...)
    return DecayChainKinematics{length(mass_tuple), length(angle_tuple), T, A}(
        SVector{length(mass_tuple), T}(mass_tuple),
        SVector{length(angle_tuple), A}(angle_tuple),
    )
end

function DecayChainKinematics(
        topology::DecayTopology,
        system::CascadeSystem;
        internal_masses2,
        vertex_angles,
    )
    masses = line_masses2(topology, system, internal_masses2)
    length(vertex_angles) == nvertices(topology) ||
        throw(ArgumentError("vertex_angles must have one entry per topology vertex"))
    return _decay_chain_kinematics_from_values(masses, vertex_angles)
end

function _require_kinematic_line_ind(x::DecayChainKinematics, line_ind::Integer)
    line_ind in Base.OneTo(length(x.line_masses2)) ||
        throw(ArgumentError("line_ind $line_ind is outside 1:$(length(x.line_masses2))"))
    return nothing
end

function _require_kinematic_vertex(x::DecayChainKinematics, vertex_ind::Integer)
    vertex_ind in Base.OneTo(length(x.vertex_angles)) ||
        throw(ArgumentError("vertex_ind $vertex_ind is outside 1:$(length(x.vertex_angles))"))
    return nothing
end

function _check_kinematics(topology::DecayTopology, x::DecayChainKinematics{Nl, Nv}) where {Nl, Nv}
    Nl == nlines(topology) ||
        throw(ArgumentError("kinematics line_masses2 must have one entry per topology line"))
    Nv == nvertices(topology) ||
        throw(ArgumentError("kinematics vertex_angles must have one entry per topology vertex"))
    return nothing
end

"""
    line_invariant(x, line_ind)
    line_invariant(topology, x, address)

Return the invariant mass squared for a topology line. Pass a line id directly,
or pass a [`DecayTopology`](@ref) and bracket address.
"""
function line_invariant(x::DecayChainKinematics, line_ind::Integer)
    _require_kinematic_line_ind(x, line_ind)
    return x.line_masses2[line_ind]
end

line_invariant(topology::DecayTopology, x::DecayChainKinematics, address) =
    line_invariant(x, line_ind_for(topology, address))

line_invariant(x::DecayChainKinematics, topology::DecayTopology, address) =
    line_invariant(topology, x, address)

"""
    vertex_masses2(topology, x, vertex_ind)
    vertex_masses2(topology, x, address)

Return `(m0², m1², m2²)` for the parent and ordered children at a topology
vertex. Pass a vertex id directly, or use a bracket address.
"""
function vertex_masses2(topology::DecayTopology, x::DecayChainKinematics, vertex_ind::Integer)
    _check_kinematics(topology, x)
    l0, l1, l2 = vertex_line_inds(topology, vertex_ind)
    return (x.line_masses2[l0], x.line_masses2[l1], x.line_masses2[l2])
end

vertex_masses2(topology::DecayTopology, x::DecayChainKinematics, address) =
    vertex_masses2(topology, x, vertex_ind_for(topology, address))

vertex_masses2(x::DecayChainKinematics, topology::DecayTopology, vertex_ind::Integer) =
    vertex_masses2(topology, x, vertex_ind)

vertex_masses2(x::DecayChainKinematics, topology::DecayTopology, address) =
    vertex_masses2(topology, x, address)

function vertex_helicities(topology::DecayTopology, two_λs, vertex_ind::Integer)
    l0, l1, l2 = vertex_line_inds(topology, vertex_ind)
    length(two_λs) == nlines(topology) ||
        throw(ArgumentError("two_λs must have one entry per topology line"))
    return (two_λs[l0], two_λs[l1], two_λs[l2])
end

function vertex_spins(topology::DecayTopology, two_js, vertex_ind::Integer)
    length(two_js) == nlines(topology) ||
        throw(ArgumentError("two_js must have one entry per topology line"))
    l0, l1, l2 = vertex_line_inds(topology, vertex_ind)
    return (two_js[l0], two_js[l1], two_js[l2])
end

"""
    vertex_angles(x, vertex_ind)
    vertex_angles(topology, x, address)

Return the local helicity-angle record `(cosθ, ϕ)` for a topology vertex. Pass a
vertex id directly, or pass a [`DecayTopology`](@ref) and bracket address.
"""
function vertex_angles(x::DecayChainKinematics, vertex_ind::Integer)
    _require_kinematic_vertex(x, vertex_ind)
    return x.vertex_angles[vertex_ind]
end

vertex_angles(topology::DecayTopology, x::DecayChainKinematics, address) =
    vertex_angles(x, vertex_ind_for(topology, address))

vertex_angles(x::DecayChainKinematics, topology::DecayTopology, address) =
    vertex_angles(topology, x, address)
