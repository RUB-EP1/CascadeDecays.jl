"""
    DecayChain(topology, spins; propagators, vertices)

Static cascade model made from a validated flat topology and typed payloads.
The public constructor accepts bracket-addressed `address => payload` pairs and
stores the resolved flat graph information in typed static arrays. The `spins`
input is consumed at construction time to build `line_two_js(chain)`; it is not
stored on the chain.
"""
struct DecayChain{
        Nf,
        Np,
        Nv,
        Nl,
        P,
        V,
        T <: DecayTopology,
    }
    topology::T
    propagators::SVector{Np, P}
    vertices::SVector{Nv, V}
    propagating_line_inds::SVector{Np, Int}
    line_two_js::SVector{Nl, Int}
end

function DecayChain(
        topology::DecayTopology,
        propagators::SVector{Np, P},
        vertices::SVector{Nv, V},
        propagating_line_inds::SVector{Np, Int},
        line_two_js::SVector{Nl, Int},
    ) where {Np, Nv, Nl, P, V}
    Nv == nvertices(topology) ||
        throw(ArgumentError("number of vertices must match topology"))
    Nl == nlines(topology) ||
        throw(ArgumentError("line_two_js must have one entry per topology line"))
    all(vertex_ind -> is_binary_vertex(topology, vertex_ind), Base.OneTo(nvertices(topology))) ||
        throw(ArgumentError("DecayChain currently requires binary topology vertices"))
    all(line_ind -> isinternal_line_ind(topology, line_ind), propagating_line_inds) ||
        throw(ArgumentError("propagators may only be attached to internal lines"))
    length(unique(propagating_line_inds)) == Np ||
        throw(ArgumentError("propagating lines must be unique"))
    sort(collect(propagating_line_inds)) == internal_line_inds(topology) ||
        throw(ArgumentError("provide exactly one propagator spec for each internal line"))
    return DecayChain{nfinal(topology), Np, Nv, Nl, P, V, typeof(topology)}(
        topology,
        propagators,
        vertices,
        propagating_line_inds,
        line_two_js,
    )
end

function DecayChain(
        topology::DecayTopology,
        spins::SystemSpinsOrSpinParities,
        propagators,
        vertices,
        propagating_line_inds,
        propagator_two_js,
    )
    propagator_tuple = Tuple(propagators)
    vertex_tuple = Tuple(vertices)
    line_tuple = Tuple(Int(line_ind) for line_ind in propagating_line_inds)
    two_j_tuple = Tuple(Int(two_j) for two_j in propagator_two_js)
    length(propagator_tuple) == length(line_tuple) ||
        throw(ArgumentError("propagators and propagating_line_inds must have the same length"))
    length(two_j_tuple) == length(line_tuple) ||
        throw(ArgumentError("propagator_two_js and propagating_line_inds must have the same length"))
    length(vertex_tuple) == nvertices(topology) ||
        throw(ArgumentError("number of vertices must match topology"))
    all(line_ind -> isinternal_line_ind(topology, line_ind), line_tuple) ||
        throw(ArgumentError("propagators may only be attached to internal lines"))
    length(unique(line_tuple)) == length(line_tuple) ||
        throw(ArgumentError("propagating lines must be unique"))
    sort(collect(line_tuple)) == internal_line_inds(topology) ||
        throw(ArgumentError("provide exactly one propagator spec for each internal line"))
    _check_spins(topology, spins)

    line_two_j_values = MVector{nlines(topology), Int}(undef)
    for (i, line_ind) in pairs(final_line_inds(topology))
        line_two_j_values[line_ind] = Int(final_two_js(spins)[i])
    end
    for (line_ind, two_j) in zip(line_tuple, two_j_tuple)
        line_two_j_values[line_ind] = Int(two_j)
    end
    line_two_j_values[root_line_ind(topology)] = Int(root_two_j(spins))

    return DecayChain(
        topology,
        SVector{length(propagator_tuple)}(propagator_tuple),
        SVector{length(vertex_tuple)}(vertex_tuple),
        SVector{length(line_tuple), Int}(line_tuple),
        SVector(line_two_j_values),
    )
end

function _vertex_payload_for(vertex_specs::Tuple, vertex_ids::Tuple, vertex_ind::Integer)
    matches = findall(==(vertex_ind), vertex_ids)
    length(matches) == 1 ||
        throw(ArgumentError("provide exactly one vertex payload for vertex_ind $vertex_ind"))
    return vertex_specs[only(matches)].second
end

"""
    DecayChain(topology, spins; propagators, vertices)

Build a static cascade model from bracket-addressed payload pairs. User-facing
addresses are resolved immediately to internal line and vertex ids. External,
internal, and root spins are stored as the line-indexed `line_two_js(chain)`;
parities are consumed by LS helpers before chain construction.
"""
function DecayChain(
        topology::DecayTopology,
        spins::SystemSpinsOrSpinParities;
        propagators::Tuple{Vararg{PropagatorSpec}},
        vertices::Tuple{Vararg{Pair{<:Any, <:Any}}},
    )
    length(vertices) == nvertices(topology) ||
        throw(ArgumentError("vertices must contain one payload per topology vertex"))
    propagating_line_tuple = Tuple(line_ind_for(topology, spec.first) for spec in propagators)
    all(line_ind -> isinternal_line_ind(topology, line_ind), propagating_line_tuple) ||
        throw(ArgumentError("propagator addresses must refer to internal lines"))
    vertex_id_tuple = Tuple(vertex_ind_for(topology, spec.first) for spec in vertices)
    length(unique(vertex_id_tuple)) == nvertices(topology) ||
        throw(ArgumentError("vertices must address each topology vertex exactly once"))

    ordered_vertices = ntuple(
        vertex_ind -> _vertex_payload_for(vertices, vertex_id_tuple, vertex_ind),
        nvertices(topology),
    )
    return DecayChain(
        topology,
        spins,
        Tuple(spec.second.lineshape for spec in propagators),
        ordered_vertices,
        propagating_line_tuple,
        Tuple(spec.second.two_j for spec in propagators),
    )
end

relation(chain::DecayChain) = relation(chain.topology)
root_line_ind(chain::DecayChain) = root_line_ind(chain.topology)
final_line_inds(chain::DecayChain) = final_line_inds(chain.topology)
nlines(chain::DecayChain) = nlines(chain.topology)
nvertices(chain::DecayChain) = nvertices(chain.topology)
nfinal(::DecayChain{Nf}) where {Nf} = Nf
internal_line_inds(chain::DecayChain) = internal_line_inds(chain.topology)

"""
    propagating_line_inds(chain)

Return the internal line ids that carry propagators, in propagator-index order.
"""
propagating_line_inds(chain::DecayChain) = chain.propagating_line_inds
bracket(chain::DecayChain; labels = nothing) = bracket(chain.topology; labels)
bracket_notation(chain::DecayChain; labels = nothing) =
    bracket_notation(chain.topology; labels)
vertex_line_inds(chain::DecayChain, vertex_ind::Integer) = vertex_line_inds(chain.topology, vertex_ind)
child_line_inds(chain::DecayChain, vertex_ind::Integer) = child_line_inds(chain.topology, vertex_ind)
final_descendants(chain::DecayChain, line_ind::Integer) = final_descendants(chain.topology, line_ind)
vertex_masses2(chain::DecayChain, x::DecayChainKinematics, vertex_ind::Integer) =
    vertex_masses2(chain.topology, x, vertex_ind)
vertex_helicities(chain::DecayChain, two_λs, vertex_ind::Integer) =
    vertex_helicities(chain.topology, two_λs, vertex_ind)

line_two_js(chain::DecayChain) = chain.line_two_js
propagator_two_js(chain::DecayChain) = chain.line_two_js[propagating_line_inds(chain)]

function vertex_spins(chain::DecayChain, vertex_ind::Integer)
    two_js = line_two_js(chain)
    l0, l1, l2 = vertex_line_inds(chain, vertex_ind)
    return (two_js[l0], two_js[l1], two_js[l2])
end
