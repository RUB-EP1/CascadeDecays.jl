module CascadeDecays

using StaticArrays
using FourVectors
using InstructionalDecayTrees
using Tullio

import PartialWaveFunctions: wignerD_doublearg
import ThreeBodyDecays
import ThreeBodyDecays: SpinParity, RecouplingLS, ⊗

# interfaces
export AbstractLineshape, AbstractVertex
export ConstantLineshape
include("interfaces.jl")

# propagators
export Propagator
include("propagator.jl")

# topology
export DecayTopology
export relation, root_line_ind, final_line_inds, nlines, nvertices, nfinal
export line_ind_for, vertex_ind_for
export incoming_line_inds, outgoing_line_inds, incoming_line_ind
export child_line_inds, vertex_line_inds, final_descendants
export produced_by, consumed_by
export internal_line_inds, propagating_line_inds
export isroot_line_ind, isfinal_line_ind, isinternal_line_ind
export has_canonical_line_order
export bracket_notation
include("topology.jl")

# quantum numbers
export SystemMasses, SystemParities, SystemSpinParities
export SystemSpins, SystemHelicities, SystemSpinProjections
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
export SpinParity, Vertex
export line_two_js
include("chain.jl")

# kinematic frames and path comparisons
export HelicityRootFrame, CurrentFrame
include("frames_to_be_moved_to_idt.jl")
export helicity_frame_path, relative_wigner_angles
include("wigner_rotations.jl")

# evaluation
export routed_vertex_amplitude, routed_propagator_product, amplitude
include("evaluation.jl")

# kinematic task / point
export KinematicTask, KinematicPoint
export kinematics_at, alignment_angles_at
include("spec_routine.jl")

# cascade container
export CascadeDecay
export reference_topology, cascade_system, couplings, names
export unpolarized_intensity
include("cascade_decay.jl")

# LS coupling
export possible_vertex_ls, minimal_vertex_coupling
export possible_vertex_couplings, minimal_vertex_couplings
export line_spin_parities, vertex_spin_parities
export minimal_ls_decay_chain, all_ls_decay_chains
include("ls_coupling.jl")

end
