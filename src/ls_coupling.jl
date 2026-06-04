using ThreeBodyDecays: SpinParity, VertexFunction, RecouplingLS, possible_ls
using StaticArrays

"""
    line_spin_parities(topology, system, propagator_specs)

Assemble a line-indexed `SpinParity` view. Requires a [`CascadeSystem`](@ref)
built from [`SystemSpinParities`](@ref). Final and root entries combine external
spins and parities; internal entries come from each [`PropagatorFunction`](@ref) built with
`PropagatorFunction(jp, lineshape)`.
"""
function line_spin_parities(
    topology::DecayTopology,
    system::CascadeSystemWithParities,
    propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
)
    _check_system(topology, system)
    parities = system.quantum.parities
    jps = MVector{nlines(topology),SpinParity}(undef)
    for (i, line) in pairs(finallines(topology))
        jps[line] = SpinParity(final_two_js(system)[i], parities.finals[i])
    end
    for spec in propagator_specs
        line = line_for(topology, spec.first)
        prop = spec.second
        jps[line] = SpinParity(prop.two_j, prop.extra.parity)
    end
    jps[rootline(topology)] = SpinParity(root_two_j(system), parities.P0)
    return SVector(jps)
end

function vertex_spin_parities(
    topology::DecayTopology,
    system::CascadeSystemWithParities,
    propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
    vertex::Integer,
)
    jps = line_spin_parities(topology, system, propagator_specs)
    l0, l1, l2 = vertex_lines(topology, vertex)
    return (jps[l0], jps[l1], jps[l2])
end

"""
    possible_vertex_ls(jp0, jp1, jp2)

Return allowed `(two_l, two_s)` couplings at one binary vertex using full J^P
constraints from `ThreeBodyDecays.possible_ls`.
"""
possible_vertex_ls(jp0::SpinParity, jp1::SpinParity, jp2::SpinParity) =
    possible_ls(jp1, jp2; jp = jp0)

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
    return ntuple(nvertices(topology)) do vertex
        vertex_address(topology, vertex) => spec_fn(vertex)
    end
end

function _vertex_coupling_values(specs)
    return map(pair -> pair.second, specs)
end

"""
    possible_vertex_couplings(topology, system, propagator_specs)

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
    system::CascadeSystemWithParities,
    propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
)
    return _vertex_coupling_specs(topology, vertex -> begin
        jp0, jp1, jp2 = vertex_spin_parities(topology, system, propagator_specs, vertex)
        possible_vertex_ls(jp0, jp1, jp2)
    end)
end

"""
    minimal_vertex_couplings(topology, system, propagator_specs)

Return a tuple of `address => (two_l, two_s)` pairs, one per topology vertex.
Pair ordering follows the same root-first preorder convention as
[`possible_vertex_couplings`](@ref).
"""
function minimal_vertex_couplings(
    topology::DecayTopology,
    system::CascadeSystemWithParities,
    propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
)
    return _vertex_coupling_specs(topology, vertex -> begin
        jp0, jp1, jp2 = vertex_spin_parities(topology, system, propagator_specs, vertex)
        minimal_vertex_coupling(jp0, jp1, jp2)
    end)
end

function _ls_vertices(vertex_couplings)
    return SVector(map(c -> VertexFunction(RecouplingLS(c)), vertex_couplings))
end

function _build_ls_decay_chain(
    topology::DecayTopology,
    propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
    vertex_couplings,
)
    propagating_line_tuple = Tuple(line_for(topology, spec.first) for spec in propagator_specs)
    return DecayChain(
        topology,
        Tuple(spec.second.lineshape for spec in propagator_specs),
        _ls_vertices(vertex_couplings),
        propagating_line_tuple,
        Tuple(spec.second.two_j for spec in propagator_specs),
    )
end

"""
    minimal_ls_decay_chain(topology, system, propagator_specs)

Build one [`DecayChain`](@ref) using the smallest allowed orbital coupling at
each vertex. `system` must be a [`CascadeSystem`](@ref) built from
[`SystemSpinParities`](@ref).

This is the same as [`all_ls_decay_chains`](@ref), but takes `first` of each
local coupling list from [`possible_vertex_couplings`](@ref).
"""
function minimal_ls_decay_chain(
    topology::DecayTopology,
    system::CascadeSystemWithParities,
    propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
)
    couplings = _vertex_coupling_values(minimal_vertex_couplings(topology, system, propagator_specs))
    return _build_ls_decay_chain(topology, propagator_specs, couplings)
end

"""
    all_ls_decay_chains(topology, system, propagator_specs)

Build every [`DecayChain`](@ref) obtained from the Cartesian product of local
LS couplings at all vertices. `system` must be a [`CascadeSystem`](@ref) built
from [`SystemSpinParities`](@ref).

What the function does is equivalent to the example below.

# Example (three-body decay)

````julia
let
    (l1, itr_v1), (l2, itr_v2) = possible_vertex_couplings(topology, weak, propagators)
    map(Iterators.product(itr_v1, itr_v2)) do ((two_l1, two_s1), (two_l2, two_s2))
        DecayChain(topology;
            propagators,
            vertices = (
                l1 => VertexFunction(RecouplingLS((two_l1, two_s1))),
                l2 => VertexFunction(RecouplingLS((two_l2, two_s2))),
            ),
        )
    end
end
````
"""
function all_ls_decay_chains(
    topology::DecayTopology,
    system::CascadeSystemWithParities,
    propagator_specs::Tuple{Vararg{PropagatorSpecWithParity}},
)
    coupling_specs = possible_vertex_couplings(topology, system, propagator_specs)
    coupling_lists = ntuple(i -> coupling_specs[i].second, nvertices(topology))
    return [
        _build_ls_decay_chain(topology, propagator_specs, combo) for
        combo in Iterators.product(coupling_lists...)
    ]
end
