module CascadeDecays

using StaticArrays

export AbstractLineshape, AbstractVertex
export ConstantLineshape
export CascadeSystem, CascadeKinematics
export DecayTopology, DecayChain
export relation, rootline, finallines, nlines, nvertices, nfinal
export incoming_lines, outgoing_lines, incoming_line
export child_lines, vertex_lines, final_descendants
export produced_by, consumed_by
export internal_lines, propagating_lines
export isrootline, isfinalline, isinternalline
export has_canonical_line_order, line_masses2
export vertex_masses2, vertex_helicities, vertex_spins, vertex_angles
export line_invariant
export bracket, validate_topology

include("interfaces.jl")
include("topology.jl")
include("system.jl")
include("kinematics.jl")
include("chain.jl")

end
