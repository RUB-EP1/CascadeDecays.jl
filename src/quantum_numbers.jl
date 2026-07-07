"""
    UndefinedParity

Parity character `'±'` when parity is not fixed or not specified. `possible_ls_more`
skips the parity filter when any participating line carries it.
"""
const UndefinedParity = Char('±')

function _spin_label(two_j::Integer)
    return iseven(two_j) ? string(two_j ÷ 2) : string(two_j, "/2")
end

function Base.show(io::IO, ::MIME"text/plain", jp::SpinParity)
    print(io, _spin_label(jp.two_j), jp.p)
    return
end

"""
    SystemMasses(m1, m2, ...; m0)

External system mass descriptor. Final-state masses are positional; the root
mass is always supplied as the keyword `m0`.
"""
struct SystemMasses{Nf, T <: Real}
    finals::SVector{Nf, T}
    m0::T
end

function SystemMasses(ms...; m0)
    mass_tuple = Tuple(ms)
    length(mass_tuple) >= 1 ||
        throw(ArgumentError("provide at least one final-state mass before `m0`"))
    T = promote_type(typeof(m0), map(typeof, mass_tuple)...)
    return SystemMasses{length(mass_tuple), T}(
        SVector{length(mass_tuple), T}(mass_tuple),
        convert(T, m0),
    )
end

function SystemMasses(ms::ThreeBodyDecays.MassTuple)
    return SystemMasses(ms.m1, ms.m2, ms.m3; m0 = ms.m0)
end

Base.:(==)(a::SystemMasses, b::SystemMasses) =
    a.m0 == b.m0 && a.finals == b.finals

function _check_helicity(two_j::Integer, two_λ::Integer)
    return abs(two_λ) <= two_j ||
        throw(ArgumentError("helicity two_λ=$two_λ is not allowed for spin two_j=$two_j"))
end

"""
    SystemSpins(two_h1, two_h2, ...; h0=nothing, two_h0=nothing)
    SystemSpins(spins, two_λ1, two_λ2, ...; h0=nothing, two_h0=nothing)

External line-indexed doubled spin or helicity descriptor. Final-state values are
positional (line `1`, `2`, …); the root is always supplied with `two_h0=...` or
`h0=...`.

Type aliases [`SystemHelicities`](@ref) and [`SystemSpinProjections`](@ref) name the
same struct for amplitude helicities and spin projections.

Pass a [`SystemSpins`](@ref) as the first argument when building helicities to check
each value against the corresponding line spin (`|two_λ| ≤ two_j`).
"""
struct SystemSpins{Nf}
    finals::SVector{Nf, Int}
    two_h0::Int
end

function _system_spins_root(; h0::Union{Nothing, Int} = nothing, two_h0::Union{Nothing, Int} = nothing)
    !isnothing(two_h0) && !isnothing(h0) &&
        throw(ArgumentError("provide root with either `two_h0=...` or `h0=...`, not both"))
    if !isnothing(two_h0)
        return two_h0
    elseif !isnothing(h0)
        return 2h0
    else
        throw(ArgumentError("provide root with `two_h0=...` or `h0=...`"))
    end
end

function SystemSpins(finals...; h0::Union{Nothing, Int} = nothing, two_h0::Union{Nothing, Int} = nothing)
    spin_tuple = Tuple(finals)
    length(spin_tuple) >= 1 ||
        throw(ArgumentError("provide at least one final-state entry before the root entry"))
    parent_two_h = _system_spins_root(; h0, two_h0)
    return SystemSpins{length(spin_tuple)}(
        SVector{length(spin_tuple), Int}(map(Int, spin_tuple)),
        parent_two_h,
    )
end

function SystemSpins(spins::SystemSpins, λs...; h0::Union{Nothing, Int} = nothing, two_h0::Union{Nothing, Int} = nothing)
    λ_tuple = Tuple(map(Int, λs))
    Nf = length(spins.finals)
    length(λ_tuple) == Nf ||
        throw(ArgumentError("provide one final-state helicity per final line in `spins`"))
    parent_two_λ = _system_spins_root(; h0, two_h0)
    _check_helicity(spins.two_h0, parent_two_λ)
    for i in 1:Nf
        _check_helicity(spins.finals[i], λ_tuple[i])
    end
    return SystemSpins(λ_tuple...; two_h0 = parent_two_λ)
end

function SystemSpins(two_js::ThreeBodyDecays.SpinTuple)
    return SystemSpins(two_js.two_h1, two_js.two_h2, two_js.two_h3; two_h0 = two_js.two_h0)
end

"""Alias of [`SystemSpins`](@ref) for external helicity assignments in [`amplitude`](@ref)."""
const SystemHelicities = SystemSpins

"""Alias of [`SystemSpins`](@ref) for external spin-projection assignments."""
const SystemSpinProjections = SystemSpins

"""
    SystemParities(P1, P2, ...; P0)

External system parity descriptor. Final-state parities are positional; the
root parity is always supplied as the keyword `P0`.
"""
struct SystemParities{Nf}
    finals::SVector{Nf, Char}
    P0::Char
end

function SystemParities(Ps...; P0)
    parity_tuple = Tuple(Ps)
    length(parity_tuple) >= 1 ||
        throw(ArgumentError("provide at least one final-state parity before `P0`"))
    return SystemParities{length(parity_tuple)}(SVector{length(parity_tuple), Char}(parity_tuple), P0)
end

function SystemParities(Ps::ThreeBodyDecays.ParityTuple)
    return SystemParities(Ps.P1, Ps.P2, Ps.P3; P0 = Ps.P0)
end

"""
    SystemSpinParities(spins, parities)
    SystemSpinParities(spins, P1, P2, ...; P0)
    SystemSpinParities(jp1, jp2, ...; jp0)
    SystemSpinParities("1+", "0-", ...; jp0="1±")

External spin and parity descriptor for [`CascadeSystem`](@ref) when LS
enumeration needs explicit final-state and root parities.

Spin-parity labels use the same `jp"…"` / string syntax as
[`ThreeBodyDecays.ThreeBodySpinParities`](https://github.com/gw2ssi/ThreeBodyDecays.jl).
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

function SystemSpinParities(
        jps::ThreeBodyDecays.SpinParity...;
        jp0::ThreeBodyDecays.SpinParity = error(
            "provide root spin-parity with `jp0`, e.g. `jp0=jp\"1±\"` or `jp0=\"1±\"`",
        ),
    )
    jp_tuple = Tuple(jps)
    length(jp_tuple) >= 1 ||
        throw(ArgumentError("provide at least one final-state jp before `jp0`"))
    spins = SystemSpins(map(jp -> jp.two_j, jp_tuple)...; two_h0 = jp0.two_j)
    parities = SystemParities(map(jp -> jp.p, jp_tuple)...; P0 = jp0.p)
    return SystemSpinParities(spins, parities)
end

function SystemSpinParities(
        jps::AbstractString...;
        jp0::AbstractString = error(
            "provide root spin-parity with `jp0`, e.g. `jp0=jp\"1±\"` or `jp0=\"1±\"`",
        ),
    )
    return SystemSpinParities(
        map(ThreeBodyDecays.str2jp, jps)...;
        jp0 = ThreeBodyDecays.str2jp(jp0),
    )
end

SystemSpins(quantum::SystemSpinParities, λs...; kwargs...) = SystemSpins(quantum.spins, λs...; kwargs...)

final_two_js(spins::SystemSpins) = spins.finals
final_two_js(quantum::SystemSpinParities) = final_two_js(quantum.spins)
root_two_j(spins::SystemSpins) = spins.two_h0
root_two_j(quantum::SystemSpinParities) = root_two_j(quantum.spins)
