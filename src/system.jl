"""
    CascadeSystem(two_js, masses2)

Static external information for a cascade model. `two_js` and `masses2` are
given as `(final_1, ..., final_n, root)`. Internal spins belong to the
`DecayChain` propagator specs.
"""
struct CascadeSystem{Nf,T<:Real}
    final_two_js::SVector{Nf,Int}
    final_masses2::SVector{Nf,T}
    root_two_j::Int
    root_mass2::T
end

function CascadeSystem(two_js, masses2)
    two_j_tuple = Tuple(Int(two_j) for two_j in two_js)
    mass_tuple = Tuple(masses2)
    length(two_j_tuple) == length(mass_tuple) ||
        throw(ArgumentError("two_js and masses2 must have the same length"))
    length(two_j_tuple) >= 2 ||
        throw(ArgumentError("two_js and masses2 must contain at least one final state and one root"))
    final_two_js = two_j_tuple[1:(end - 1)]
    final_masses2 = mass_tuple[1:(end - 1)]
    root_two_j = two_j_tuple[end]
    root_mass2 = mass_tuple[end]
    T = promote_type(typeof(root_mass2), map(typeof, final_masses2)...)
    return CascadeSystem{length(final_masses2),T}(
        SVector{length(final_two_js),Int}(final_two_js),
        SVector{length(final_masses2),T}(final_masses2),
        Int(root_two_j),
        convert(T, root_mass2),
    )
end

function _check_system(topology::DecayTopology, system::CascadeSystem{Nf}) where {Nf}
    Nf == nfinal(topology) ||
        throw(ArgumentError("system final_masses2 must have one entry per final line"))
    return nothing
end

"""
    line_values(topology; finals, root, internals=())

Assemble a complete line-indexed `SVector` from public topology addresses.
`finals` are ordered by final-state labels, `root` is the mother value, and
`internals` are bracket-addressed pairs such as `(1, 2) => value`.
"""
function line_values(topology::DecayTopology; finals, root, internals = ())
    final_tuple = Tuple(finals)
    internal_specs = Tuple(internals)
    length(final_tuple) == nfinal(topology) ||
        throw(ArgumentError("finals must have one entry per final-state line"))
    all(spec -> spec isa Pair, internal_specs) ||
        throw(ArgumentError("internals must be provided as `address => value` pairs"))
    internal_line_tuple = Tuple(line_for(topology, spec.first) for spec in internal_specs)
    all(line -> isinternalline(topology, line), internal_line_tuple) ||
        throw(ArgumentError("internal addresses must refer to internal lines"))
    sort(collect(internal_line_tuple)) == internal_lines(topology) ||
        throw(ArgumentError("provide exactly one internal value for each internal line"))

    T = promote_type(typeof(root), map(typeof, final_tuple)..., map(spec -> typeof(spec.second), internal_specs)...)
    values = MVector{nlines(topology),T}(undef)
    for (i, line) in pairs(finallines(topology))
        values[line] = final_tuple[i]
    end
    for (spec, line) in zip(internal_specs, internal_line_tuple)
        values[line] = spec.second
    end
    values[rootline(topology)] = root
    return SVector(values)
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
