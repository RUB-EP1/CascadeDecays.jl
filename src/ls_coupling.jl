"""
    line_spin_parities(topology, quantum, propagator_specs)

Assemble a line-indexed `SpinParity` view. Requires [`SystemSpinParities`](@ref).
Final and root entries combine external spins and parities; internal entries come
from each [`Propagator`](@ref) built with `Propagator(jp, lineshape)`.
"""
function line_spin_parities(
        topology::DecayTopology,
        quantum::SystemSpinParities,
        propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
    )
    _check_spins(topology, quantum.spins)
    parities = quantum.parities
    jps = MVector{nlines(topology), SpinParity}(undef)
    for (i, line_ind) in pairs(final_line_inds(topology))
        jps[line_ind] = SpinParity(final_two_js(quantum)[i], parities.finals[i])
    end
    for spec in propagator_specs
        line_ind = line_ind_for(topology, spec.first)
        prop = spec.second
        jps[line_ind] = SpinParity(prop.two_j, prop.extra.parity)
    end
    jps[root_line_ind(topology)] = SpinParity(root_two_j(quantum), parities.P0)
    return SVector(jps)
end

function vertex_spin_parities(
        topology::DecayTopology,
        quantum::SystemSpinParities,
        propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
        vertex_ind::Integer,
    )
    jps = line_spin_parities(topology, quantum, propagator_specs)
    l0, l1, l2 = vertex_line_inds(topology, vertex_ind)
    return (jps[l0], jps[l1], jps[l2])
end

"""
    possible_vertex_ls(jp0, jp1, jp2)

Return allowed `(two_l, two_s)` couplings at one binary vertex using
`possible_ls_more` (including [`UndefinedParity`](@ref), `'±'`).
"""
possible_vertex_ls(jp0::SpinParity, jp1::SpinParity, jp2::SpinParity) =
    possible_ls_more(jp1, jp2; jp = jp0)

"""
    minimal_vertex_coupling(jp0, jp1, jp2)

Smallest-`two_l` coupling allowed by [`possible_vertex_ls`](@ref).
"""
minimal_vertex_coupling(jp0, jp1, jp2) = _first_coupling(possible_vertex_ls(jp0, jp1, jp2))

function _first_coupling(couplings)
    isempty(couplings) &&
        throw(ArgumentError("no LS couplings allowed for the supplied J^P assignment"))
    return first(couplings)
end

function _vertex_coupling_specs(topology::DecayTopology, spec_fn)
    return ntuple(nvertices(topology)) do vertex_ind
        vertex_address(topology, vertex_ind) => spec_fn(vertex_ind)
    end
end

function _vertex_coupling_values(specs)
    return map(pair -> pair.second, specs)
end

"""
    possible_vertex_couplings(topology, quantum, propagator_specs)

Return a tuple of `address => couplings` pairs, one per topology vertex.
Each `address` is a bracket key such as `(1, 2)` or `((1, 2), 3)`; each
`couplings` is a list of allowed `(two_l, two_s)` tuples sorted by increasing
`two_l`.

Pairs appear in **root-first preorder**: the outer/root-attached vertex first,
then deeper vertices from left to right — the same order as
[`DecayChain`](@ref) expects for its `vertices` keyword. Within each list,
couplings are sorted by increasing `two_l`.
"""
function possible_vertex_couplings(
        topology::DecayTopology,
        quantum::SystemSpinParities,
        propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
    )
    return _vertex_coupling_specs(
        topology, vertex_ind -> begin
            jp0, jp1, jp2 = vertex_spin_parities(topology, quantum, propagator_specs, vertex_ind)
            possible_vertex_ls(jp0, jp1, jp2)
        end
    )
end

"""
    minimal_vertex_couplings(topology, quantum, propagator_specs)

Return a tuple of `address => (two_l, two_s)` pairs, one per topology vertex.
Pair ordering follows the same root-first preorder convention as
[`possible_vertex_couplings`](@ref).
"""
function minimal_vertex_couplings(
        topology::DecayTopology,
        quantum::SystemSpinParities,
        propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
    )
    return _vertex_coupling_specs(
        topology, vertex_ind -> begin
            jp0, jp1, jp2 = vertex_spin_parities(topology, quantum, propagator_specs, vertex_ind)
            minimal_vertex_coupling(jp0, jp1, jp2)
        end
    )
end

function _ls_vertices(vertex_couplings)
    return SVector(map(c -> Vertex(RecouplingLS(c)), vertex_couplings))
end

function _build_ls_decay_chain(
        topology::DecayTopology,
        quantum::SystemSpinParities,
        propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
        vertex_couplings,
    )
    propagating_line_tuple = Tuple(line_ind_for(topology, spec.first) for spec in propagator_specs)
    return DecayChain(
        topology,
        quantum.spins,
        Tuple(spec.second.lineshape for spec in propagator_specs),
        _ls_vertices(vertex_couplings),
        propagating_line_tuple,
        Tuple(spec.second.two_j for spec in propagator_specs),
    )
end

"""
    minimal_ls_decay_chain(topology, quantum, propagator_specs)

Build one [`DecayChain`](@ref) using the smallest allowed orbital coupling at
each vertex. `quantum` must be a [`SystemSpinParities`](@ref) value.

This is the same as [`all_ls_decay_chains`](@ref), but takes `first` of each
local coupling list from [`possible_vertex_couplings`](@ref).
"""
function minimal_ls_decay_chain(
        topology::DecayTopology,
        quantum::SystemSpinParities,
        propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
    )
    couplings = _vertex_coupling_values(minimal_vertex_couplings(topology, quantum, propagator_specs))
    return _build_ls_decay_chain(topology, quantum, propagator_specs, couplings)
end

"""
    all_ls_decay_chains(topology, quantum, propagator_specs)

Build every [`DecayChain`](@ref) obtained from the Cartesian product of local
LS couplings at all vertices. `quantum` must be a [`SystemSpinParities`](@ref)
value.

What the function does is equivalent to the example below.

# Example (three-body decay)

````julia
let
    (l1, itr_v1), (l2, itr_v2) = possible_vertex_couplings(topology, weak, propagators)
    map(Iterators.product(itr_v1, itr_v2)) do ((two_l1, two_s1), (two_l2, two_s2))
        DecayChain(topology, weak.spins;
            propagators,
            vertices = (
                l1 => Vertex(RecouplingLS((two_l1, two_s1))),
                l2 => Vertex(RecouplingLS((two_l2, two_s2))),
            ),
        )
    end
end
````
"""
function all_ls_decay_chains(
        topology::DecayTopology,
        quantum::SystemSpinParities,
        propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
    )
    coupling_specs = possible_vertex_couplings(topology, quantum, propagator_specs)
    coupling_lists = ntuple(i -> coupling_specs[i].second, nvertices(topology))
    return [
        _build_ls_decay_chain(topology, quantum, propagator_specs, combo) for
            combo in Iterators.product(coupling_lists...)
    ]
end
