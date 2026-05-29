using ThreeBodyDecays: SpinParity, VertexFunction, RecouplingLS, possible_ls
using StaticArrays

"""
    ParityAugmentedSystem(system, parities)

Spin system from [`CascadeSystem`](@ref) paired with explicit [`SystemParities`](@ref).
Construct with [`add_parities`](@ref) for automatic LS chain builders.
"""
struct ParityAugmentedSystem{Nf,Tm}
    system::CascadeSystem{Nf,Tm}
    parities::SystemParities{Nf}
end

"""
    add_parities(system, parities)
    add_parities(system, P1, P2, ...; P0)

Attach external and root parities to a [`CascadeSystem`](@ref) built from
[`SystemSpins`](@ref). Use [`UndefinedParity`](@ref) on `P0` when the mother
parity is unknown (typical weak decays); final-state parities stay explicit.

The result is accepted by [`minimal_ls_decay_chain`](@ref) and
[`all_ls_decay_chains`](@ref).
"""
function add_parities(system::CascadeSystem{Nf,Tm}, parities::SystemParities{Nf}) where {Nf,Tm}
    length(parities.finals) == Nf ||
        throw(ArgumentError("parities must have one entry per final-state line"))
    return ParityAugmentedSystem(system, parities)
end

add_parities(system::CascadeSystem{Nf,Tm}, Ps...; P0) where {Nf,Tm} =
    add_parities(system, SystemParities(Ps...; P0))

"""
    line_spin_parities(topology, system_with_parities, propagator_specs)

Assemble a line-indexed `SpinParity` view. Final and root entries combine spins
from the underlying [`CascadeSystem`](@ref) with [`SystemParities`](@ref);
internal entries come from each propagator spec's `jp` (or `two_j` with `p`).
"""
function line_spin_parities(topology::DecayTopology, system::ParityAugmentedSystem, propagator_specs)
    _check_system(topology, system.system)
    spec_tuple = _require_pair_specs(:propagators, propagator_specs)
    jps = MVector{nlines(topology),SpinParity}(undef)
    for (i, line) in pairs(finallines(topology))
        jps[line] = SpinParity(final_two_js(system.system)[i], system.parities.finals[i])
    end
    for spec in spec_tuple
        line = line_for(topology, spec.first)
        jps[line] = _propagator_spin_parity(spec.second)
    end
    jps[rootline(topology)] = SpinParity(root_two_j(system.system), system.parities.P0)
    return SVector(jps)
end

function vertex_spin_parities(topology::DecayTopology, system::ParityAugmentedSystem, propagator_specs, vertex::Integer)
    jps = line_spin_parities(topology, system, propagator_specs)
    l0, l1, l2 = vertex_lines(topology, vertex)
    return (jps[l0], jps[l1], jps[l2])
end

vertex_spin_parities(topology::DecayTopology, system::ParityAugmentedSystem, propagator_specs, address) =
    vertex_spin_parities(topology, system, propagator_specs, vertex_for(topology, address))

_spin_parity_defined(jp::SpinParity) = parity_defined(jp.p)

function _possible_ls_spin_only(two_j0::Integer, two_j1::Integer, two_j2::Integer)
    two_ls = Tuple{Int,Int}[]
    for two_s in abs(two_j1 - two_j2):2:(two_j1 + two_j2)
        for two_l in abs(two_j0 - two_s):2:(two_j0 + two_s)
            push!(two_ls, (two_l, two_s))
        end
    end
    return sort(two_ls, by = x -> x[1])
end

"""
    possible_vertex_ls(jp0, jp1, jp2)

Return allowed `(two_l, two_s)` couplings at one binary vertex. When any entry
has [`UndefinedParity`](@ref), parity constraints are skipped at that vertex.
"""
function possible_vertex_ls(jp0::SpinParity, jp1::SpinParity, jp2::SpinParity)
    if _spin_parity_defined(jp0) && _spin_parity_defined(jp1) && _spin_parity_defined(jp2)
        return possible_ls(jp1, jp2; jp = jp0)
    end
    return _possible_ls_spin_only(jp0.two_j, jp1.two_j, jp2.two_j)
end

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
    possible_vertex_couplings(topology, system_with_parities, propagator_specs)

Return a tuple of `address => couplings` pairs, one per topology vertex.
Each `address` is a bracket key such as `(1, 2)` or `((1, 2), 3)`; each
`couplings` is a list of allowed `(two_l, two_s)` tuples sorted by increasing
`two_l`.

Pairs appear in **root-first preorder**: the outer/root-attached vertex first,
then deeper vertices from left to right — the same order as
[`DecayChain`](@ref) expects for its `vertices` keyword. Within each list,
couplings are sorted by increasing `two_l`.
"""
function possible_vertex_couplings(topology::DecayTopology, system::ParityAugmentedSystem, propagator_specs)
    return _vertex_coupling_specs(topology, vertex -> begin
        jp0, jp1, jp2 = vertex_spin_parities(topology, system, propagator_specs, vertex)
        possible_vertex_ls(jp0, jp1, jp2)
    end)
end

"""
    minimal_vertex_couplings(topology, system_with_parities, propagator_specs)

Return a tuple of `address => (two_l, two_s)` pairs, one per topology vertex.
Pair ordering follows the same root-first preorder convention as
[`possible_vertex_couplings`](@ref).
"""
function minimal_vertex_couplings(topology::DecayTopology, system::ParityAugmentedSystem, propagator_specs)
    return _vertex_coupling_specs(topology, vertex -> begin
        jp0, jp1, jp2 = vertex_spin_parities(topology, system, propagator_specs, vertex)
        minimal_vertex_coupling(jp0, jp1, jp2)
    end)
end

function _ls_vertices(vertex_couplings)
    return SVector(map(c -> VertexFunction(RecouplingLS(c)), vertex_couplings))
end

function _build_ls_decay_chain(topology::DecayTopology, propagator_specs, vertex_couplings)
    spec_tuple = _require_pair_specs(:propagators, propagator_specs)
    propagating_line_tuple = Tuple(line_for(topology, spec.first) for spec in spec_tuple)
    return DecayChain(
        topology,
        Tuple(_propagator_lineshape(spec.second) for spec in spec_tuple),
        _ls_vertices(vertex_couplings),
        propagating_line_tuple,
        Tuple(_propagator_two_j(spec.second) for spec in spec_tuple),
    )
end

"""
    minimal_ls_decay_chain(topology, system_with_parities, propagator_specs)

Build one [`DecayChain`](@ref) using the smallest allowed orbital coupling at
each vertex. Pass [`add_parities`](@ref)(`system`, parities) as the second
argument.

This is the same as [`all_ls_decay_chains`](@ref), but takes `first` of each
local coupling list from [`possible_vertex_couplings`](@ref).
"""
function minimal_ls_decay_chain(topology::DecayTopology, system::ParityAugmentedSystem, propagator_specs)
    couplings = _vertex_coupling_values(minimal_vertex_couplings(topology, system, propagator_specs))
    return _build_ls_decay_chain(topology, propagator_specs, couplings)
end

"""
    all_ls_decay_chains(topology, system_with_parities, propagator_specs)

Build every [`DecayChain`](@ref) obtained from the Cartesian product of local
LS couplings at all vertices. Pass [`add_parities`](@ref)(`system`, parities)
as the second argument.

What the function does is equivalent to the example below.

# Example (three-body decay)

```julia
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
```
"""
function all_ls_decay_chains(topology::DecayTopology, system::ParityAugmentedSystem, propagator_specs)
    coupling_specs = possible_vertex_couplings(topology, system, propagator_specs)
    coupling_lists = ntuple(i -> coupling_specs[i].second, nvertices(topology))
    return [
        _build_ls_decay_chain(topology, propagator_specs, combo) for
        combo in Iterators.product(coupling_lists...)
    ]
end
