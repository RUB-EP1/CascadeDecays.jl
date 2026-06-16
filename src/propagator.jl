"""
    Propagator(two_j, lineshape)
    Propagator(jp, lineshape)

Propagator payload for an internal cascade line, analogous to
[`Vertex`](@ref) for vertices.

- `two_j`: twice the spin of the propagating system
- `lineshape`: callable propagator model (typically [`AbstractLineshape`](@ref))
- `extra`: optional metadata (`nothing`, or `(; parity=...)` from `SpinParity`)

See also [`PropagatorWithParity`](@ref) for the LS-coupling case.
"""
struct Propagator{F, X}
    two_j::Int
    lineshape::F
    extra::X
end

const PropagatorParityExtra = NamedTuple{(:parity,)}

"""Alias of [`Propagator`](@ref) with parity metadata in `extra`."""
const PropagatorWithParity{F} = Propagator{F, <:PropagatorParityExtra}
const PropagatorSpec = Pair{<:Tuple, <:Propagator}
const PropagatorSpecWithParity = Pair{<:Tuple, <:PropagatorWithParity}

function Propagator(two_j::Integer, lineshape::F) where {F}
    return Propagator{F, Nothing}(Int(two_j), lineshape, nothing)
end

function Propagator(jp::SpinParity, lineshape::F) where {F}
    return Propagator{F, NamedTuple{(:parity,), Tuple{typeof(jp.p)}}}(
        jp.two_j,
        lineshape,
        (parity = jp.p,),
    )
end
