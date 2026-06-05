"""
    AbstractLineshape

Abstract supertype for propagator or lineshape payloads attached to internal
cascade lines.
"""
abstract type AbstractLineshape end

"""
    AbstractVertex

Abstract supertype for local binary-decay vertex payloads.
"""
abstract type AbstractVertex end

"""
    Vertex{R,F}

Local two-body decay vertex payload (alias of `ThreeBodyDecays.VertexFunction`).

Wraps a recoupling scheme `h::R` and optional form factor `ff::F`. Used in
[`DecayChain`](@ref) vertex slots.

See also [`AbstractVertex`](@ref).
"""
const Vertex = ThreeBodyDecays.VertexFunction

"""
    ConstantLineshape(value)

Minimal callable lineshape useful for tests and early prototypes.
"""
struct ConstantLineshape{T} <: AbstractLineshape
    value::T
end

(p::ConstantLineshape)(σ) = p.value
