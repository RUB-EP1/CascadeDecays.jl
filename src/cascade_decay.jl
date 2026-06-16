"""
    CascadeDecay(chains, system, reference_topology; couplings)

Container for a coherent set of decay chains that share one system and reference
topology.
"""
struct CascadeDecay{Nc,S,T,C}
    chains::NTuple{Nc,DecayChain}
    couplings::SVector{Nc,C}
    system::S
    reference_topology::T
end

function CascadeDecay(
    chains::Tuple{Vararg{DecayChain}},
    system::CascadeSystem,
    reference_topology::DecayTopology,
    couplings::Tuple{Vararg{Number}},
)
    isempty(chains) && throw(ArgumentError("chains must contain at least one DecayChain"))
    length(couplings) == length(chains) ||
        throw(ArgumentError("couplings must have one entry per chain"))
    coupling_tuple = Tuple(couplings)
    C = promote_type(map(typeof, coupling_tuple)...)
    return CascadeDecay{length(chains),typeof(system),typeof(reference_topology),C}(
        chains,
        SVector{length(coupling_tuple),C}(coupling_tuple),
        system,
        reference_topology,
    )
end

function CascadeDecay(
    chains::Tuple{Vararg{DecayChain}},
    system::CascadeSystem,
    reference_topology::DecayTopology;
    couplings::Tuple{Vararg{Number}}=ntuple(_ -> one(ComplexF64), length(chains)),
)
    return CascadeDecay(chains, system, reference_topology, couplings)
end

reference_topology(cascade::CascadeDecay) = cascade.reference_topology
cascade_system(cascade::CascadeDecay) = cascade.system
couplings(cascade::CascadeDecay) = cascade.couplings
