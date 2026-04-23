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
    ConstantLineshape(value)

Minimal callable lineshape useful for tests and early prototypes.
"""
struct ConstantLineshape{T} <: AbstractLineshape
    value::T
end

(p::ConstantLineshape)(σ) = p.value
