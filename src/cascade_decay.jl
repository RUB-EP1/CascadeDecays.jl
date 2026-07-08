"""
    CascadeDecay(chains, reference_topology; couplings, names)

Container for a coherent set of decay chains that share compatible external
spin axes and one reference topology. Each chain has a coupling coefficient and
a string name for selection and display.

`chains` is stored as a concretely typed tuple of [`DecayChain`](@ref) objects.
`couplings` and `names` are homogeneous `NTuple`s of one numeric type and
`String`, respectively. Vector inputs are accepted at the constructor boundary
and converted immediately.
"""
struct CascadeDecay{Nc, Ch <: Tuple{Vararg{DecayChain}}, T, C}
    chains::Ch
    couplings::NTuple{Nc, C}
    names::NTuple{Nc, String}
    reference_topology::T
end

function _external_two_js(chain::DecayChain)
    two_js = line_two_js(chain)
    return (Tuple(two_js[final_line_inds(chain)])..., two_js[root_line_ind(chain)])
end

function _validate_cascade_inputs(
        reference_topology::DecayTopology,
        chains::Tuple{Vararg{DecayChain}},
    )
    expected_nfinal = nfinal(reference_topology)
    expected_external_two_js = _external_two_js(chains[1])
    length(expected_external_two_js) == expected_nfinal + 1 ||
        throw(ArgumentError("chains must have the same number of final particles as reference_topology"))
    for chain in chains
        nfinal(chain) == expected_nfinal ||
            throw(ArgumentError("all chains must have the same number of final particles as reference_topology"))
        _external_two_js(chain) == expected_external_two_js ||
            throw(ArgumentError("all chains in a CascadeDecay must share external spin assignments"))
    end
    return nothing
end

function _default_chain_name(_chain::DecayChain, index::Integer)
    return "chain_$index"
end

function _default_chain_names(chains::Tuple{Vararg{DecayChain}})
    return ntuple(i -> _default_chain_name(chains[i], i), length(chains))
end

function _validate_chain_names(chains::Tuple, names::Tuple)
    return length(names) == length(chains) ||
        throw(ArgumentError("names must have one entry per chain"))
end

function _coupling_tuple(couplings::Tuple{Vararg{Number}})
    C = promote_type(map(typeof, couplings)...)
    return ntuple(i -> convert(C, couplings[i]), length(couplings))
end

function _name_tuple(names::Tuple{Vararg{AbstractString}})
    return ntuple(i -> String(names[i]), length(names))
end

function CascadeDecay(
        chains::Tuple{Vararg{DecayChain}},
        reference_topology::DecayTopology,
        couplings::Tuple{Vararg{Number}},
        names::Tuple{Vararg{AbstractString}},
    )
    isempty(chains) && throw(ArgumentError("chains must contain at least one DecayChain"))
    length(couplings) == length(chains) ||
        throw(ArgumentError("couplings must have one entry per chain"))
    _validate_chain_names(chains, names)
    _validate_cascade_inputs(reference_topology, chains)
    coupling_tuple = _coupling_tuple(couplings)
    names_tuple = _name_tuple(names)
    C = eltype(coupling_tuple)
    return CascadeDecay{length(chains), typeof(chains), typeof(reference_topology), C}(
        chains,
        coupling_tuple,
        names_tuple,
        reference_topology,
    )
end

function CascadeDecay(
        chains::Tuple{Vararg{DecayChain}},
        reference_topology::DecayTopology,
        couplings::Tuple{Vararg{Number}};
        names::Union{Nothing, Tuple{Vararg{AbstractString}}} = nothing,
    )
    name_tuple = isnothing(names) ? _default_chain_names(chains) : names
    return CascadeDecay(chains, reference_topology, couplings, name_tuple)
end

function CascadeDecay(
        chains::Tuple{Vararg{DecayChain}},
        reference_topology::DecayTopology;
        couplings::Tuple{Vararg{Number}} = ntuple(_ -> one(ComplexF64), length(chains)),
        names::Union{Nothing, Tuple{Vararg{AbstractString}}} = nothing,
    )
    return CascadeDecay(chains, reference_topology, couplings; names)
end

function CascadeDecay(
        chains::AbstractVector{<:DecayChain},
        reference_topology::DecayTopology,
        couplings::AbstractVector{<:Number},
        names::AbstractVector{<:AbstractString},
    )
    return CascadeDecay(
        (chains...,),
        reference_topology,
        (couplings...,),
        (String(name) for name in names...),
    )
end

function CascadeDecay(
        chains::AbstractVector{<:DecayChain},
        reference_topology::DecayTopology,
        couplings::AbstractVector{<:Number};
        names::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    )
    name_tuple = isnothing(names) ? _default_chain_names((chains...,)) : (String(name) for name in names...)
    return CascadeDecay((chains...,), reference_topology, (couplings...,), name_tuple)
end

function CascadeDecay(
        chains::AbstractVector{<:DecayChain},
        reference_topology::DecayTopology;
        couplings::AbstractVector{<:Number} = fill(one(ComplexF64), length(chains)),
        names::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    )
    return CascadeDecay(chains, reference_topology, couplings; names)
end

reference_topology(cascade::CascadeDecay) = cascade.reference_topology
couplings(cascade::CascadeDecay) = cascade.couplings

Base.length(cascade::CascadeDecay) = length(cascade.chains)

function Base.iterate(cascade::CascadeDecay, state = 1)
    state > length(cascade) && return nothing
    return (
        (
            chain = cascade.chains[state],
            coupling = cascade.couplings[state],
            name = cascade.names[state],
        ),
        state + 1,
    )
end

function _slice_cascade(cascade::CascadeDecay, inds::AbstractVector{<:Integer})
    isempty(inds) && throw(ArgumentError("cannot slice CascadeDecay to zero chains"))
    n = length(inds)
    chain_tuple = ntuple(i -> cascade.chains[inds[i]], n)
    coupling_tuple = ntuple(i -> cascade.couplings[inds[i]], n)
    names_tuple = ntuple(i -> cascade.names[inds[i]], n)
    return CascadeDecay(
        chain_tuple,
        cascade.reference_topology,
        coupling_tuple,
        names_tuple,
    )
end

Base.getindex(cascade::CascadeDecay, inds::AbstractVector{Bool}) =
    _slice_cascade(cascade, findall(inds))

Base.getindex(cascade::CascadeDecay, inds::Tuple{Vararg{Bool}}) =
    _slice_cascade(cascade, findall(identity, collect(inds)))

Base.getindex(cascade::CascadeDecay, inds::AbstractVector{<:Integer}) =
    _slice_cascade(cascade, collect(inds))

Base.getindex(cascade::CascadeDecay, inds::UnitRange{<:Integer}) =
    _slice_cascade(cascade, collect(inds))

function Base.getindex(cascade::CascadeDecay, inds::NTuple{N, Integer}) where {N}
    chain_tuple = ntuple(i -> cascade.chains[inds[i]], Val(N))
    coupling_tuple = ntuple(i -> cascade.couplings[inds[i]], Val(N))
    names_tuple = ntuple(i -> cascade.names[inds[i]], Val(N))
    return CascadeDecay(
        chain_tuple,
        cascade.reference_topology,
        coupling_tuple,
        names_tuple,
    )
end

function Base.getindex(cascade::CascadeDecay, ind::Integer)
    return CascadeDecay(
        (cascade.chains[ind],),
        cascade.reference_topology,
        (cascade.couplings[ind],),
        (cascade.names[ind],),
    )
end

function Base.getindex(cascade::CascadeDecay, name::AbstractString)
    inds = findall(==(name), cascade.names)
    isempty(inds) && throw(KeyError(name))
    return cascade[inds]
end

function Base.getindex(cascade::CascadeDecay, names::AbstractVector{<:AbstractString})
    inds = Int[]
    for name in names
        matches = findall(==(name), cascade.names)
        isempty(matches) && throw(KeyError(name))
        append!(inds, matches)
    end
    return cascade[inds]
end

"""
    merge(cascade1, cascade2, cascades...)

Concatenate decay chains from compatible [`CascadeDecay`](@ref) models. All
inputs must share the same `reference_topology` and compatible external spin
axes. Duplicate chain names are rejected.

This extends `Base.merge` for `CascadeDecay` arguments.
"""
function Base.merge(cascade1::CascadeDecay, cascade2::CascadeDecay, cascades::CascadeDecay...)
    return foldl(merge, (cascade2, cascades...); init = cascade1)
end

function Base.merge(cascade1::CascadeDecay, cascade2::CascadeDecay)
    cascade1.reference_topology == cascade2.reference_topology ||
        throw(ArgumentError("merged models must share reference_topology"))
    _external_two_js(cascade1.chains[1]) == _external_two_js(cascade2.chains[1]) ||
        throw(ArgumentError("merged models must share external spin assignments"))
    combined_names = (cascade1.names..., cascade2.names...)
    length(combined_names) == length(unique(combined_names)) ||
        throw(ArgumentError("merged models contain duplicate chain names"))
    chains = (cascade1.chains..., cascade2.chains...)
    couplings = (cascade1.couplings..., cascade2.couplings...)
    return CascadeDecay(
        chains,
        cascade1.reference_topology,
        couplings,
        combined_names,
    )
end

function _chain_topology_string(chain::DecayChain)
    return bracket_notation(chain.topology)
end

function _coupling_string(coupling)
    coupling isa Complex && return string(round(real(coupling); digits = 6), " + ", round(imag(coupling); digits = 6), "im")
    return sprint(show, coupling)
end

function Base.show(io::IO, cascade::CascadeDecay)
    print(io, "CascadeDecay with ", length(cascade), " chain", length(cascade) == 1 ? "" : "s")
    print(io, " on ", bracket_notation(cascade.reference_topology), ":")
    name_width = max(4, maximum(length(name) for name in cascade.names; init = 0))
    coupling_width = max(
        8,
        maximum(length(_coupling_string(c)) for c in cascade.couplings; init = 0),
    )
    topology_width = max(
        8,
        maximum(length(_chain_topology_string(row.chain)) for row in cascade; init = 0),
    )
    print(io, "\n  ", lpad("name", name_width), "  ", lpad("coupling", coupling_width), "  ", lpad("topology", topology_width))
    for row in cascade
        print(
            io,
            "\n  ",
            lpad(row.name, name_width),
            "  ",
            lpad(_coupling_string(row.coupling), coupling_width),
            "  ",
            lpad(_chain_topology_string(row.chain), topology_width),
        )
    end
    return
end

"""
    amplitude(cascade, point::KinematicPoint)

Coherent helicity amplitude `sum(cᵢ * Aᵢ)` for all chains in `cascade`.
"""
function amplitude(cascade::CascadeDecay{Nc}, point::KinematicPoint) where {Nc}
    point.task.reference_topology == cascade.reference_topology ||
        throw(ArgumentError("point task reference_topology must match cascade reference_topology"))
    amp1 = amplitude(cascade.chains[1], point)
    Nc == 1 && return cascade.couplings[1] .* amp1
    res = cascade.couplings[1] .* amp1
    for i in 2:Nc
        res .+= cascade.couplings[i] .* amplitude(cascade.chains[i], point)
    end
    return res
end

"""
    unpolarized_intensity(cascade, point::KinematicPoint)

Return `sum(abs2, amplitude(cascade, point))`.
"""
unpolarized_intensity(cascade::CascadeDecay, point::KinematicPoint) =
    sum(abs2, amplitude(cascade, point))
