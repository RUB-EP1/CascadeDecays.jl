"""
    PropagatorFunction(two_j, lineshape)
    PropagatorFunction(jp, lineshape)

Propagator payload for an internal cascade line, analogous to
`ThreeBodyDecays.VertexFunction` for vertices.

- `two_j`: twice the spin of the propagating system
- `lineshape`: callable propagator model (typically [`AbstractLineshape`](@ref))
- `extra`: optional metadata (`nothing`, or `(; parity=...)` from `SpinParity`)

See also [`PropagatorFunctionWithParity`](@ref) for the LS-coupling case.
"""
struct PropagatorFunction{F,X}
    two_j::Int
    lineshape::F
    extra::X
end

const PropagatorParityExtra = NamedTuple{(:parity,)}

"""Alias of [`PropagatorFunction`](@ref) with parity metadata in `extra`."""
const PropagatorFunctionWithParity{F} = PropagatorFunction{F,<:PropagatorParityExtra}
const PropagatorSpec = Pair{<:Tuple,<:PropagatorFunction}
const PropagatorSpecWithParity = Pair{<:Tuple,<:PropagatorFunctionWithParity}

function PropagatorFunction(two_j::Integer, lineshape::F) where {F}
    return PropagatorFunction{F,Nothing}(Int(two_j), lineshape, nothing)
end

function PropagatorFunction(jp::SpinParity, lineshape::F) where {F}
    return PropagatorFunction{F,NamedTuple{(:parity,),Tuple{typeof(jp.p)}}}(
        jp.two_j,
        lineshape,
        (parity = jp.p,),
    )
end
