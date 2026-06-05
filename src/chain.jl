"""
    DecayChain(topology; propagators, vertices)

Static cascade model made from a validated flat topology and typed payloads.
The public constructor accepts bracket-addressed `address => payload` pairs and
stores the resolved flat graph information in typed static arrays.
"""
struct DecayChain{
    Nf,
    Np,
    Nv,
    P,
    V,
    T<:DecayTopology,
}
    topology::T
    propagators::SVector{Np,P}
    vertices::SVector{Nv,V}
    propagating_lines::SVector{Np,Int}
    propagator_two_js::SVector{Np,Int}
end

function DecayChain(
    topology::DecayTopology,
    propagators::SVector{Np,P},
    vertices::SVector{Nv,V},
    propagating_lines::SVector{Np,Int},
    propagator_two_js::SVector{Np,Int},
) where {Np,Nv,P,V}
    Nv == nvertices(topology) ||
        throw(ArgumentError("number of vertices must match topology"))
    all(line -> isinternalline(topology, line), propagating_lines) ||
        throw(ArgumentError("propagators may only be attached to internal lines"))
    length(unique(propagating_lines)) == Np ||
        throw(ArgumentError("propagating lines must be unique"))
    sort(collect(propagating_lines)) == internal_lines(topology) ||
        throw(ArgumentError("provide exactly one propagator spec for each internal line"))
    return DecayChain{nfinal(topology),Np,Nv,P,V,typeof(topology)}(
        topology,
        propagators,
        vertices,
        propagating_lines,
        propagator_two_js,
    )
end

function DecayChain(topology::DecayTopology, propagators, vertices, propagating_lines, propagator_two_js)
    propagator_tuple = Tuple(propagators)
    vertex_tuple = Tuple(vertices)
    line_tuple = Tuple(Int(line) for line in propagating_lines)
    two_j_tuple = Tuple(Int(two_j) for two_j in propagator_two_js)
    length(propagator_tuple) == length(line_tuple) ||
        throw(ArgumentError("propagators and propagating_lines must have the same length"))
    length(two_j_tuple) == length(line_tuple) ||
        throw(ArgumentError("propagator_two_js and propagating_lines must have the same length"))
    length(vertex_tuple) == nvertices(topology) ||
        throw(ArgumentError("number of vertices must match topology"))
    return DecayChain(
        topology,
        SVector{length(propagator_tuple)}(propagator_tuple),
        SVector{length(vertex_tuple)}(vertex_tuple),
        SVector{length(line_tuple),Int}(line_tuple),
        SVector{length(two_j_tuple),Int}(two_j_tuple),
    )
end

function _vertex_payload_for(vertex_specs::Tuple, vertex_ids::Tuple, vertex_ind::Integer)
    matches = findall(==(vertex_ind), vertex_ids)
    length(matches) == 1 ||
        throw(ArgumentError("provide exactly one vertex payload for vertex_ind $vertex_ind"))
    return vertex_specs[only(matches)].second
end

"""
    DecayChain(topology; propagators, vertices)

Build a static cascade model from bracket-addressed payload pairs. User-facing
addresses are resolved immediately to internal line and vertex ids.
"""
function DecayChain(
    topology::DecayTopology;
    propagators::Tuple{Vararg{PropagatorSpec}},
    vertices::Tuple{Vararg{Pair{<:Any,<:Any}}},
)
    length(vertices) == nvertices(topology) ||
        throw(ArgumentError("vertices must contain one payload per topology vertex"))
    propagating_line_tuple = Tuple(line_for(topology, spec.first) for spec in propagators)
    all(line -> isinternalline(topology, line), propagating_line_tuple) ||
        throw(ArgumentError("propagator addresses must refer to internal lines"))
    vertex_id_tuple = Tuple(vertex_for(topology, spec.first) for spec in vertices)
    length(unique(vertex_id_tuple)) == nvertices(topology) ||
        throw(ArgumentError("vertices must address each topology vertex exactly once"))

    ordered_vertices = ntuple(
        vertex_ind -> _vertex_payload_for(vertices, vertex_id_tuple, vertex_ind),
        nvertices(topology),
    )
    return DecayChain(
        topology,
        Tuple(spec.second.lineshape for spec in propagators),
        ordered_vertices,
        propagating_line_tuple,
        Tuple(spec.second.two_j for spec in propagators),
    )
end

relation(chain::DecayChain) = relation(chain.topology)
rootline(chain::DecayChain) = rootline(chain.topology)
finallines(chain::DecayChain) = finallines(chain.topology)
nlines(chain::DecayChain) = nlines(chain.topology)
nvertices(chain::DecayChain) = nvertices(chain.topology)
nfinal(::DecayChain{Nf}) where {Nf} = Nf
internal_lines(chain::DecayChain) = internal_lines(chain.topology)
propagating_lines(chain::DecayChain) = chain.propagating_lines
bracket(chain::DecayChain; labels = nothing) = bracket(chain.topology; labels)
vertex_lines(chain::DecayChain, vertex_ind::Integer) = vertex_lines(chain.topology, vertex_ind)
child_lines(chain::DecayChain, vertex_ind::Integer) = child_lines(chain.topology, vertex_ind)
final_descendants(chain::DecayChain, line::Integer) = final_descendants(chain.topology, line)
vertex_masses2(chain::DecayChain, x::CascadeKinematics, vertex_ind::Integer) =
    vertex_masses2(chain.topology, x, vertex_ind)
vertex_helicities(chain::DecayChain, two_λs, vertex_ind::Integer) =
    vertex_helicities(chain.topology, two_λs, vertex_ind)

function line_two_js(chain::DecayChain, system::CascadeSystem)
    _check_system(chain.topology, system)
    spins = MVector{nlines(chain),Int}(undef)
    for (i, line) in pairs(finallines(chain))
        spins[line] = final_two_js(system)[i]
    end
    for (i, line) in pairs(propagating_lines(chain))
        spins[line] = chain.propagator_two_js[i]
    end
    spins[rootline(chain)] = root_two_j(system)
    return SVector(spins)
end

function vertex_spins(chain::DecayChain, system::CascadeSystem, vertex_ind::Integer)
    two_js = line_two_js(chain, system)
    l0, l1, l2 = vertex_lines(chain, vertex_ind)
    return (two_js[l0], two_js[l1], two_js[l2])
end
