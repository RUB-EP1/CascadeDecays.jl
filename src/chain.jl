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
    propagating_line_inds::SVector{Np,Int}
    propagator_two_js::SVector{Np,Int}
end

function DecayChain(
    topology::DecayTopology,
    propagators::SVector{Np,P},
    vertices::SVector{Nv,V},
    propagating_line_inds::SVector{Np,Int},
    propagator_two_js::SVector{Np,Int},
) where {Np,Nv,P,V}
    Nv == nvertices(topology) ||
        throw(ArgumentError("number of vertices must match topology"))
    all(line_ind -> isinternal_line_ind(topology, line_ind), propagating_line_inds) ||
        throw(ArgumentError("propagators may only be attached to internal lines"))
    length(unique(propagating_line_inds)) == Np ||
        throw(ArgumentError("propagating lines must be unique"))
    sort(collect(propagating_line_inds)) == internal_line_inds(topology) ||
        throw(ArgumentError("provide exactly one propagator spec for each internal line"))
    return DecayChain{nfinal(topology),Np,Nv,P,V,typeof(topology)}(
        topology,
        propagators,
        vertices,
        propagating_line_inds,
        propagator_two_js,
    )
end

function DecayChain(topology::DecayTopology, propagators, vertices, propagating_line_inds, propagator_two_js)
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
vertex_masses2(chain::DecayChain, x::CascadeKinematics, vertex_ind::Integer) =
    vertex_masses2(chain.topology, x, vertex_ind)
vertex_helicities(chain::DecayChain, two_λs, vertex_ind::Integer) =
    vertex_helicities(chain.topology, two_λs, vertex_ind)

function line_two_js(chain::DecayChain{Nf,Np}, system::CascadeSystem) where {Nf,Np}
    _check_system(chain.topology, system)
    final_inds = Tuple(final_line_inds(chain))
    prop_inds = Tuple(propagating_line_inds(chain))
    final_js = final_two_js(system)
    prop_js = chain.propagator_two_js
    root_j = root_two_j(system)
    return SVector(ntuple(Val(Nf + Np + 1)) do line_ind
        i = findfirst(==(line_ind), final_inds)
        i !== nothing && return Int(final_js[i])
        i = findfirst(==(line_ind), prop_inds)
        i !== nothing && return Int(prop_js[i])
        return Int(root_j)
    end)
end

function vertex_spins(chain::DecayChain, system::CascadeSystem, vertex_ind::Integer)
    two_js = line_two_js(chain, system)
    l0, l1, l2 = vertex_line_inds(chain, vertex_ind)
    return (two_js[l0], two_js[l1], two_js[l2])
end
