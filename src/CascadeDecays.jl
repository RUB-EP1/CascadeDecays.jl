module CascadeDecays

using StaticArrays

export AbstractLineshape, AbstractVertex
export ConstantLineshape
export DecayTopology, DecayChain
export relation, rootline, finallines, nlines, nvertices, nfinal
export incoming_lines, outgoing_lines, incoming_line
export produced_by, consumed_by
export internal_lines, propagating_lines
export isrootline, isfinalline, isinternalline
export bracket, validate_topology

include("interfaces.jl")
include("topology.jl")
include("chain.jl")

end
