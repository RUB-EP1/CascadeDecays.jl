module CascadeDecays

using StaticArrays
using FourVectors
using InstructionalDecayTrees

import ThreeBodyDecays
import ThreeBodyDecays: SpinParity, VertexFunction, RecouplingLS, ⊗

# interfaces
export AbstractLineshape, AbstractVertex
export ConstantLineshape
include("interfaces.jl")

# propagators
export PropagatorFunction
include("propagator_function.jl")

# topology
export DecayTopology
export relation, rootline, finallines, nlines, nvertices, nfinal
export incoming_lines, outgoing_lines, incoming_line
export child_lines, vertex_lines, final_descendants
export produced_by, consumed_by
export internal_lines, propagating_lines
export isrootline, isfinalline, isinternalline
export has_canonical_line_order
export bracket, validate_topology
include("topology.jl")

# quantum numbers
export SystemMasses, SystemSpins, SystemParities, SystemSpinParities
export UndefinedParity
export final_two_js, root_two_j
include("quantum_numbers.jl")

# possible_ls_more (internal; used by ls_coupling)
include("possible_ls_more.jl")

# system
export CascadeSystem, CascadeSystemWithParities
export final_masses, root_mass
export line_masses2, line_values
include("system.jl")

# kinematics
export CascadeKinematics
export line_invariant
export vertex_masses2, vertex_helicities, vertex_spins, vertex_angles
include("kinematics.jl")

# decay chain
export DecayChain
export SpinParity, VertexFunction
export line_two_js
include("chain.jl")

# evaluation
export HelicityRootFrame, CurrentFrame
export helicity_angle_program, helicity_angle_programs, cascade_kinematics
export routed_vertex_amplitude, routed_propagator_product, amplitude
include("evaluation.jl")

# LS coupling
export possible_vertex_ls, minimal_vertex_coupling
export possible_vertex_couplings, minimal_vertex_couplings
export line_spin_parities, vertex_spin_parities
export minimal_ls_decay_chain, all_ls_decay_chains
include("ls_coupling.jl")

end
