"""
    CascadeSystem(spins, masses)
    CascadeSystem(spin_parity, masses)

Static external information for a cascade model. Pass [`SystemSpins`](@ref) for
spin-only systems, or [`SystemSpinParity`](@ref) when LS builders need explicit
external and root parities. Masses are always [`SystemMasses`](@ref).
"""
struct CascadeSystem{Nf,Tm}
    spins::SystemSpins{Nf}
    masses::SystemMasses{Nf,Tm}
    parities::Union{SystemParities{Nf}, Nothing}
end

function CascadeSystem(spins::SystemSpins{Nf}, masses::SystemMasses{Nf,Tm}) where {Nf,Tm}
    return CascadeSystem{Nf,Tm}(spins, masses, nothing)
end

function CascadeSystem(quantum::SystemSpinParity{Nf}, masses::SystemMasses{Nf,Tm}) where {Nf,Tm}
    return CascadeSystem{Nf,Tm}(quantum.spins, masses, quantum.parities)
end

has_parities(system::CascadeSystem) = !isnothing(system.parities)

function _require_system_parities(system::CascadeSystem)
    has_parities(system) ||
        throw(ArgumentError("LS builders require CascadeSystem(SystemSpinParity(...), masses)"))
    return system.parities
end

final_two_js(system::CascadeSystem) = final_two_js(system.spins)
root_two_j(system::CascadeSystem) = root_two_j(system.spins)
final_masses(system::CascadeSystem) = system.masses.finals
root_mass(system::CascadeSystem) = system.masses.m0

function _check_system(topology::DecayTopology, system::CascadeSystem{Nf}) where {Nf}
    Nf == nfinal(topology) ||
        throw(ArgumentError("system final_masses must have one entry per final line"))
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

Assemble a complete line-indexed invariant-mass-squared view. External and
root masses come from `system` and are squared on assembly; internal entries
are supplied directly as runtime invariants in `internal_masses2`.
"""
function line_masses2(topology::DecayTopology, system::CascadeSystem, internal_masses2)
    _check_system(topology, system)
    internal_tuple = Tuple(internal_masses2)
    internal = internal_lines(topology)
    length(internal_tuple) == length(internal) ||
        throw(ArgumentError("internal_masses2 must match the number of internal lines"))

    invariant_type(m) = typeof(m * m)
    T = promote_type(
        invariant_type(root_mass(system)),
        map(invariant_type, Tuple(final_masses(system)))...,
        map(typeof, internal_tuple)...,
    )
    masses = MVector{nlines(topology),T}(undef)
    for (i, line) in pairs(finallines(topology))
        masses[line] = final_masses(system)[i]^2
    end
    for (i, line) in pairs(internal)
        masses[line] = internal_tuple[i]
    end
    masses[rootline(topology)] = root_mass(system)^2
    return SVector(masses)
end
