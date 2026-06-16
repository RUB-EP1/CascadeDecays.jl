"""
    CascadeDecay(chains, system, reference_topology; couplings, names)

Container for a coherent set of decay chains that share one system and reference
topology. Each chain has a coupling coefficient and a string name for selection
and display.
"""
struct CascadeDecay{Nc,Ch<:Tuple{Vararg{DecayChain}},S,T,C}
    chains::Ch
    couplings::SVector{Nc,C}
    names::SVector{Nc,String}
    system::S
    reference_topology::T
end

function _chain_tuple(chains::Union{Tuple{Vararg{DecayChain}},AbstractVector{<:DecayChain}})
    return chains isa Tuple ? chains : (chains...,)
end

function _coupling_tuple(couplings::Union{Tuple{Vararg{Number}},AbstractVector{<:Number}})
    return couplings isa Tuple ? couplings : (couplings...,)
end

function _name_tuple(names::Union{Tuple{Vararg{AbstractString}},AbstractVector{<:AbstractString}})
    return Tuple(String(name) for name in names)
end

function _validate_cascade_inputs(
    reference_topology::DecayTopology,
    system::CascadeSystem,
    chains::Tuple{Vararg{DecayChain}},
)
    _check_system(reference_topology, system)
    for chain in chains
        _check_system(chain.topology, system)
    end
    return nothing
end

function _default_chain_name(_chain::DecayChain, index::Integer)
    return "chain_$index"
end

function _default_chain_names(chains::Tuple{Vararg{DecayChain}})
    return ntuple(i -> _default_chain_name(chains[i], i), length(chains))
end

function _validate_chain_names(chains, names)
    length(names) == length(chains) ||
        throw(ArgumentError("names must have one entry per chain"))
end

function CascadeDecay(
    chains::Union{Tuple{Vararg{DecayChain}},AbstractVector{<:DecayChain}},
    system::CascadeSystem,
    reference_topology::DecayTopology,
    couplings::Union{Tuple{Vararg{Number}},AbstractVector{<:Number}},
    names::Union{Tuple{Vararg{AbstractString}},AbstractVector{<:AbstractString}},
)
    chain_tuple = _chain_tuple(chains)
    coupling_tuple = _coupling_tuple(couplings)
    names_tuple = _name_tuple(names)
    isempty(chain_tuple) && throw(ArgumentError("chains must contain at least one DecayChain"))
    length(coupling_tuple) == length(chain_tuple) ||
        throw(ArgumentError("couplings must have one entry per chain"))
    _validate_chain_names(chain_tuple, names_tuple)
    _validate_cascade_inputs(reference_topology, system, chain_tuple)
    C = promote_type(map(typeof, coupling_tuple)...)
    return CascadeDecay{length(chain_tuple),typeof(chain_tuple),typeof(system),typeof(reference_topology),C}(
        chain_tuple,
        SVector{length(coupling_tuple),C}(coupling_tuple),
        SVector{length(names_tuple),String}(names_tuple),
        system,
        reference_topology,
    )
end

function CascadeDecay(
    chains::Tuple{Vararg{DecayChain}},
    system::CascadeSystem,
    reference_topology::DecayTopology,
    couplings::Tuple{Vararg{Number}};
    names::Union{Nothing,Tuple{Vararg{AbstractString}}}=nothing,
)
    name_tuple = isnothing(names) ? _default_chain_names(chains) : names
    return CascadeDecay(chains, system, reference_topology, couplings, name_tuple)
end

function CascadeDecay(
    chains::AbstractVector{<:DecayChain},
    system::CascadeSystem,
    reference_topology::DecayTopology,
    couplings::AbstractVector{<:Number};
    names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
)
    name_tuple = isnothing(names) ? _default_chain_names(_chain_tuple(chains)) : names
    return CascadeDecay(chains, system, reference_topology, couplings, name_tuple)
end

function CascadeDecay(
    chains::Tuple{Vararg{DecayChain}},
    system::CascadeSystem,
    reference_topology::DecayTopology;
    couplings::Tuple{Vararg{Number}}=ntuple(_ -> one(ComplexF64), length(chains)),
    names::Union{Nothing,Tuple{Vararg{AbstractString}}}=nothing,
)
    return CascadeDecay(chains, system, reference_topology, couplings; names)
end

function CascadeDecay(
    chains::AbstractVector{<:DecayChain},
    system::CascadeSystem,
    reference_topology::DecayTopology;
    couplings::AbstractVector{<:Number}=fill(one(ComplexF64), length(chains)),
    names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
)
    return CascadeDecay(chains, system, reference_topology, couplings; names)
end

reference_topology(cascade::CascadeDecay) = cascade.reference_topology
cascade_system(cascade::CascadeDecay) = cascade.system
couplings(cascade::CascadeDecay) = cascade.couplings
names(cascade::CascadeDecay) = cascade.names

Base.length(cascade::CascadeDecay) = length(cascade.chains)

function Base.iterate(cascade::CascadeDecay, state=1)
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
    chain_tuple = Tuple(cascade.chains[i] for i in inds)
    return CascadeDecay(
        chain_tuple,
        cascade.system,
        cascade.reference_topology,
        cascade.couplings[inds],
        cascade.names[inds],
    )
end

Base.getindex(cascade::CascadeDecay, inds::AbstractVector{Bool}) =
    _slice_cascade(cascade, findall(inds))

Base.getindex(cascade::CascadeDecay, inds::AbstractVector{<:Integer}) =
    _slice_cascade(cascade, collect(inds))

Base.getindex(cascade::CascadeDecay, inds::UnitRange{<:Integer}) =
    _slice_cascade(cascade, collect(inds))

function Base.getindex(cascade::CascadeDecay, inds::NTuple{N,Integer}) where {N}
    chain_tuple = ntuple(i -> cascade.chains[inds[i]], Val(N))
    coupling_tuple = ntuple(i -> cascade.couplings[inds[i]], Val(N))
    names_tuple = ntuple(i -> cascade.names[inds[i]], Val(N))
    return CascadeDecay(
        chain_tuple,
        cascade.system,
        cascade.reference_topology,
        coupling_tuple,
        names_tuple,
    )
end

function Base.getindex(cascade::CascadeDecay, ind::Integer)
    return CascadeDecay(
        (cascade.chains[ind],),
        cascade.system,
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
inputs must share the same `system` and `reference_topology`. Duplicate chain
names are rejected.

This extends `Base.merge` for `CascadeDecay` arguments.
"""
function Base.merge(cascade1::CascadeDecay, cascade2::CascadeDecay, cascades::CascadeDecay...)
    return foldl(merge, (cascade2, cascades...); init = cascade1)
end

function Base.merge(cascade1::CascadeDecay, cascade2::CascadeDecay)
    cascade1.system == cascade2.system &&
        cascade1.reference_topology == cascade2.reference_topology ||
        throw(ArgumentError("merged models must share system and reference_topology"))
    combined_names = vcat(cascade1.names, cascade2.names)
    length(combined_names) == length(unique(combined_names)) ||
        throw(ArgumentError("merged models contain duplicate chain names"))
    chains = (cascade1.chains..., cascade2.chains...)
    couplings = (cascade1.couplings..., cascade2.couplings...)
    return CascadeDecay(
        chains,
        cascade1.system,
        cascade1.reference_topology,
        couplings,
        Tuple(combined_names),
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
end

"""
    amplitude(cascade, point::KinematicPoint)

Coherent helicity amplitude `sum(cᵢ * Aᵢ)` for all chains in `cascade`.
"""
function amplitude(cascade::CascadeDecay{Nc}, point::KinematicPoint) where {Nc}
    point.task.reference_topology == cascade.reference_topology ||
        throw(ArgumentError("point task reference_topology must match cascade reference_topology"))
    amp1 = amplitude(cascade.chains[1], cascade.system, point)
    Nc == 1 && return cascade.couplings[1] .* amp1
    res = cascade.couplings[1] .* amp1
    for i in 2:Nc
        res .+= cascade.couplings[i] .* amplitude(cascade.chains[i], cascade.system, point)
    end
    return res
end

"""
    unpolarized_intensity(cascade, point::KinematicPoint)

Return `sum(abs2, amplitude(cascade, point))`.
"""
unpolarized_intensity(cascade::CascadeDecay, point::KinematicPoint) =
    sum(abs2, amplitude(cascade, point))
