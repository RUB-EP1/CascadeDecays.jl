module CascadeDecays

using StaticArrays
using FourVectors
using InstructionalDecayTrees

import ThreeBodyDecays
import ThreeBodyDecays: SpinParity

export AbstractLineshape, AbstractVertex
export ConstantLineshape
export SystemMasses, SystemSpins, SystemParities
export UndefinedParity, parity_defined
export final_two_js, root_two_j, final_masses, root_mass
export CascadeSystem, CascadeKinematics
export ParityAugmentedSystem, add_parities
export SpinParity
export DecayTopology, DecayChain
export relation, rootline, finallines, nlines, nvertices, nfinal
export incoming_lines, outgoing_lines, incoming_line
export child_lines, vertex_lines, final_descendants
export produced_by, consumed_by
export internal_lines, propagating_lines
export isrootline, isfinalline, isinternalline
export has_canonical_line_order, line_masses2
export line_values, line_two_js
export vertex_masses2, vertex_helicities, vertex_spins, vertex_angles
export line_invariant
export HelicityRootFrame, CurrentFrame
export helicity_angle_program, helicity_angle_programs, cascade_kinematics
export routed_vertex_amplitude, routed_propagator_product, amplitude
export bracket, validate_topology
export possible_vertex_ls, minimal_vertex_coupling
export possible_vertex_couplings, minimal_vertex_couplings
export line_spin_parities, vertex_spin_parities
export minimal_ls_decay_chain, all_ls_decay_chains

include("interfaces.jl")
include("topology.jl")
include("quantum_numbers.jl")
include("system.jl")
include("kinematics.jl")
include("chain.jl")
include("evaluation.jl")
include("ls_coupling.jl")

end
