"""
    CascadeSystem(two_js, final_masses2; root_mass2)

Static line-level information for a cascade model.

The current convention is that all line-indexed views contain final-state lines,
internal lines, and the root/mother line. `two_js` is indexed by canonical line id.
`final_masses2` is indexed in the order `finallines(topology)`, while
`root_mass2` stores the fixed mother/root invariant mass squared.
"""
struct CascadeSystem{Nl,Nf,T<:Real}
    two_js::SVector{Nl,Int}
    final_masses2::SVector{Nf,T}
    root_mass2::T
end

function CascadeSystem(
    two_js,
    final_masses2;
    root_mass2,
)
    two_j_tuple = Tuple(Int(two_j) for two_j in two_js)
    final_mass_tuple = Tuple(final_masses2)
    T = promote_type(typeof(root_mass2), map(typeof, final_mass_tuple)...)
    return CascadeSystem{length(two_j_tuple),length(final_mass_tuple),T}(
        SVector{length(two_j_tuple),Int}(two_j_tuple),
        SVector{length(final_mass_tuple),T}(final_mass_tuple),
        convert(T, root_mass2),
    )
end

function _check_system(topology::DecayTopology, system::CascadeSystem{Nl,Nf}) where {Nl,Nf}
    Nl == nlines(topology) ||
        throw(ArgumentError("system two_js must have one entry per topology line"))
    Nf == nfinal(topology) ||
        throw(ArgumentError("system final_masses2 must have one entry per final line"))
    return nothing
end

"""
    line_masses2(topology, system, internal_masses2)

Assemble a complete line-indexed mass-squared view. Final and root masses come
from `system`; internal masses are supplied in the order `internal_lines(topology)`.
"""
function line_masses2(topology::DecayTopology, system::CascadeSystem, internal_masses2)
    _check_system(topology, system)
    internal_tuple = Tuple(internal_masses2)
    internal = internal_lines(topology)
    length(internal_tuple) == length(internal) ||
        throw(ArgumentError("internal_masses2 must match the number of internal lines"))

    T = promote_type(
        typeof(system.root_mass2),
        eltype(system.final_masses2),
        map(typeof, internal_tuple)...,
    )
    masses = MVector{nlines(topology),T}(undef)
    for (i, line) in pairs(finallines(topology))
        masses[line] = system.final_masses2[i]
    end
    for (i, line) in pairs(internal)
        masses[line] = internal_tuple[i]
    end
    masses[rootline(topology)] = system.root_mass2
    return SVector(masses)
end
