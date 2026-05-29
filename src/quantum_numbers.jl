"""
    UndefinedParity

Character used when parity is not specified. LS enumeration then ignores parity
constraints and keeps only spin-coupling rules.
"""
const UndefinedParity = '?'

parity_defined(p::Char) = p != UndefinedParity

"""
    SystemMasses(m1, m2, ...; m0)

External system mass descriptor. Final-state masses are positional; the root
mass is always supplied as the keyword `m0`.
"""
struct SystemMasses{Nf,T<:Real}
    finals::SVector{Nf,T}
    m0::T
end

function SystemMasses(ms...; m0)
    mass_tuple = Tuple(ms)
    length(mass_tuple) >= 1 ||
        throw(ArgumentError("provide at least one final-state mass before `m0`"))
    T = promote_type(typeof(m0), map(typeof, mass_tuple)...)
    return SystemMasses{length(mass_tuple),T}(
        SVector{length(mass_tuple),T}(mass_tuple),
        convert(T, m0),
    )
end

function SystemMasses(ms::ThreeBodyDecays.MassTuple)
    return SystemMasses(ms.m1, ms.m2, ms.m3; m0 = ms.m0)
end

Base.:(==)(a::SystemMasses, b::SystemMasses) =
    a.m0 == b.m0 && a.finals == b.finals

"""
    SystemSpins(two_h1, two_h2, ...; h0=nothing, two_h0=nothing)

External system spin descriptor. Final-state spins are positional; the root
spin is always supplied with `two_h0=...` or `h0=...`.
"""
struct SystemSpins{Nf}
    finals::SVector{Nf,Int}
    two_h0::Int
end

function SystemSpins(finals...; h0 = nothing, two_h0 = nothing)
    spin_tuple = Tuple(finals)
    length(spin_tuple) >= 1 ||
        throw(ArgumentError("provide at least one final-state spin before the root spin"))
    isnothing(h0) && isnothing(two_h0) &&
        throw(ArgumentError("provide root spin with `two_h0=...` or `h0=...`"))
    parent_two_h = isnothing(two_h0) ? Int(2h0) : Int(two_h0)
    return SystemSpins{length(spin_tuple)}(
        SVector{length(spin_tuple),Int}(map(Int, spin_tuple)),
        parent_two_h,
    )
end

function SystemSpins(two_js::ThreeBodyDecays.SpinTuple)
    return SystemSpins(two_js.two_h1, two_js.two_h2, two_js.two_h3; two_h0 = two_js.two_h0)
end

"""
    SystemParities(P1, P2, ...; P0)

External system parity descriptor. Final-state parities are positional; the
root parity is always supplied as the keyword `P0`.
"""
struct SystemParities{Nf}
    finals::SVector{Nf,Char}
    P0::Char
end

function SystemParities(Ps...; P0)
    parity_tuple = Tuple(Ps)
    length(parity_tuple) >= 1 ||
        throw(ArgumentError("provide at least one final-state parity before `P0`"))
    return SystemParities{length(parity_tuple)}(SVector{length(parity_tuple),Char}(parity_tuple), P0)
end

function SystemParities(Ps::ThreeBodyDecays.ParityTuple)
    return SystemParities(Ps.P1, Ps.P2, Ps.P3; P0 = Ps.P0)
end

"""
    SystemSpinParities(spins, parities)
    SystemSpinParities(spins, P1, P2, ...; P0)

External spin and parity descriptor for [`CascadeSystem`](@ref) when LS
enumeration needs explicit final-state and root parities.
"""
struct SystemSpinParities{Nf}
    spins::SystemSpins{Nf}
    parities::SystemParities{Nf}

    function SystemSpinParities(spins::SystemSpins{Nf}, parities::SystemParities{Nf}) where {Nf}
        length(parities.finals) == Nf ||
            throw(ArgumentError("parities must have one entry per final-state line"))
        return new{Nf}(spins, parities)
    end
end

SystemSpinParities(spins::SystemSpins, Ps...; P0) =
    SystemSpinParities(spins, SystemParities(Ps...; P0))

final_two_js(spins::SystemSpins) = spins.finals
final_two_js(quantum::SystemSpinParities) = final_two_js(quantum.spins)
root_two_j(spins::SystemSpins) = spins.two_h0
root_two_j(quantum::SystemSpinParities) = root_two_j(quantum.spins)
