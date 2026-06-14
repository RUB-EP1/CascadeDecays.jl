"""
    CascadeDecay(chains, system, reference_topology; couplings)

Coherent multi-chain container: amplitudes are summed as
`sum(cᵢ * amplitude(chainᵢ, ...))`.
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

"""
    amplitude(cascade, point::KinematicPoint)

Coherent helicity amplitude `sum(cᵢ * Aᵢ)` aligned to the task reference topology.
"""
function amplitude(cascade::CascadeDecay, point::KinematicPoint)
    if !isempty(point.task.wigner_finals) &&
       point.task.reference_topology != cascade.reference_topology
        throw(ArgumentError("task reference_topology must match cascade.reference_topology"))
    end
    return sum(
        c * amplitude(chain, cascade.system, point)
        for (c, chain) in zip(cascade.couplings, cascade.chains)
    )
end

"""
    amplitude(cascade, kinematics)

Coherent helicity amplitude without reference-topology alignment.
"""
function amplitude(
    cascade::CascadeDecay,
    kinematics::Tuple{Vararg{CascadeKinematics}},
)
    length(kinematics) == length(cascade.chains) ||
        throw(ArgumentError("kinematics must have one entry per chain"))
    return sum(
        c * amplitude(chain, cascade.system, x)
        for (c, chain, x) in zip(cascade.couplings, cascade.chains, kinematics)
    )
end

function amplitude(
    cascade::CascadeDecay,
    chain_ind::Integer,
    x::CascadeKinematics,
)
    chain = cascade.chains[chain_ind]
    return amplitude(chain, cascade.system, x)
end

function amplitude(
    cascade::CascadeDecay,
    kinematics::Tuple{Vararg{CascadeKinematics}},
    external_two_λs::SystemSpins,
)
    return sum(
        c * amplitude(chain, cascade.system, x, external_two_λs)
        for (c, chain, x) in zip(cascade.couplings, cascade.chains, kinematics)
    )
end

unpolarized_intensity(cascade::CascadeDecay, kinematics::Tuple{Vararg{CascadeKinematics}}) =
    sum(abs2, amplitude(cascade, kinematics))

unpolarized_intensity(cascade::CascadeDecay, point::KinematicPoint) =
    sum(abs2, amplitude(cascade, point))
