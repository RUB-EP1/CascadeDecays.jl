abstract type AbstractInitialFrame end

"""
    HelicityRootFrame()

Start generated helicity-angle programs by transforming to the root/system
helicity frame. This is the default for fully general four-vectors.
"""
struct HelicityRootFrame <: AbstractInitialFrame end

"""
    CurrentFrame()

Start generated helicity-angle programs in the current axes. Use this when the
input four-vectors are already expressed in the intended system frame.
"""
struct CurrentFrame <: AbstractInitialFrame end

_initial_frame_program(topology::DecayTopology, ::HelicityRootFrame) =
    (ToHelicityFrame(_indices_for_line_ind(topology, root_line_ind(topology))),)
_initial_frame_program(::DecayTopology, ::CurrentFrame) = ()
