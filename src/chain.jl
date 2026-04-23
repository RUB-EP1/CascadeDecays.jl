"""
    DecayChain(topology, propagators, vertices, propagating_lines)

Static cascade model made from a validated flat topology and typed payloads.
Propagators are attached only to internal lines listed by `propagating_lines`.
"""
struct DecayChain{
    Nf,
    Np,
    Nv,
    P<:AbstractLineshape,
    V<:AbstractVertex,
    T<:DecayTopology,
}
    topology::T
    propagators::SVector{Np,P}
    vertices::SVector{Nv,V}
    propagating_lines::SVector{Np,Int}
end

function DecayChain(
    topology::DecayTopology,
    propagators::SVector{Np,P},
    vertices::SVector{Nv,V},
    propagating_lines::SVector{Np,Int},
) where {Np,Nv,P<:AbstractLineshape,V<:AbstractVertex}
    Nv == nvertices(topology) ||
        throw(ArgumentError("number of vertices must match topology"))
    all(line -> isinternalline(topology, line), propagating_lines) ||
        throw(ArgumentError("propagators may only be attached to internal lines"))
    length(unique(propagating_lines)) == Np ||
        throw(ArgumentError("propagating lines must be unique"))
    return DecayChain{nfinal(topology),Np,Nv,P,V,typeof(topology)}(
        topology,
        propagators,
        vertices,
        propagating_lines,
    )
end

function DecayChain(topology::DecayTopology, propagators, vertices, propagating_lines)
    propagator_tuple = Tuple(propagators)
    vertex_tuple = Tuple(vertices)
    line_tuple = Tuple(Int(line) for line in propagating_lines)
    return DecayChain(
        topology,
        SVector{length(propagator_tuple)}(propagator_tuple),
        SVector{length(vertex_tuple)}(vertex_tuple),
        SVector{length(line_tuple),Int}(line_tuple),
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
