function _check_spins(topology::DecayTopology, spins::SystemSpinsOrSpinParities)
    length(final_two_js(spins)) == nfinal(topology) ||
        throw(ArgumentError("spins must have one entry per final line"))
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
    internal_line_tuple = Tuple(line_ind_for(topology, spec.first) for spec in internal_specs)
    all(line_ind -> isinternal_line_ind(topology, line_ind), internal_line_tuple) ||
        throw(ArgumentError("internal addresses must refer to internal lines"))
    sort(collect(internal_line_tuple)) == internal_line_inds(topology) ||
        throw(ArgumentError("provide exactly one internal value for each internal line"))

    T = promote_type(typeof(root), map(typeof, final_tuple)..., map(spec -> typeof(spec.second), internal_specs)...)
    values = MVector{nlines(topology), T}(undef)
    for (i, line_ind) in pairs(final_line_inds(topology))
        values[line_ind] = final_tuple[i]
    end
    for (spec, line_ind) in zip(internal_specs, internal_line_tuple)
        values[line_ind] = spec.second
    end
    values[root_line_ind(topology)] = root
    return SVector(values)
end
