"""
    DecayTopology(bracket)

Construct a decay topology from bracket notation, for example
`DecayTopology((((1, 2), 3), 4))` or `DecayTopology((1, 2, 3, 4))`.

Final-state particles are numbered by the integer leaves. The left/right order
inside each tuple is preserved because helicity conventions depend on child
order. Internally the topology is stored as a flat line-vertex graph.
"""
struct DecayTopology{Nl, Nv, Nf, T <: Integer, L, C}
    relation::SMatrix{Nl, Nv, T, L}
    root::Int
    finals::SVector{Nf, Int}
    child_order::C
end

function _decay_topology_from_relation(
        relation::SMatrix{Nl, Nv, T, L},
        root::Integer,
        finals::SVector{Nf, Int};
        child_order::C,
        validate::Bool = true,
    ) where {Nl, Nv, Nf, T <: Integer, L, C}
    topology = DecayTopology{Nl, Nv, Nf, T, L, C}(relation, Int(root), finals, child_order)
    validate && validate_topology(topology)
    return topology
end

function _decay_topology_from_relation(
        relation::AbstractMatrix{T};
        root::Integer,
        finals,
        child_order,
        validate::Bool = true,
    ) where {T <: Integer}
    Nl, Nv = size(relation)
    final_tuple = Tuple(Int(line_ind) for line_ind in finals)
    final_lines = SVector{length(final_tuple), Int}(final_tuple)
    static_relation = SMatrix{Nl, Nv, T, Nl * Nv}(relation)
    static_child_order = _normalize_child_order(child_order, Nv)
    return _decay_topology_from_relation(static_relation, root, final_lines; child_order = static_child_order, validate)
end

function _normalize_child_order(child_order, nv::Integer)
    if child_order isa AbstractMatrix
        size(child_order, 2) == nv ||
            throw(ArgumentError("child_order matrix must have one column per vertex"))
        N = size(child_order, 1)
        return ntuple(v -> ntuple(row -> Int(child_order[row, v]), Val(N)), nv)
    end
    tuple_order = Tuple(child_order)
    length(tuple_order) == nv ||
        throw(ArgumentError("child_order must have one entry per vertex"))
    if !isempty(tuple_order) && all(x -> length(x) == length(first(tuple_order)), tuple_order)
        N = length(first(tuple_order))
        return ntuple(v -> NTuple{N, Int}(tuple_order[v]), nv)
    end
    return ntuple(v -> Tuple(Int(line_ind) for line_ind in tuple_order[v]), nv)
end

_bracket_leaf(x::Integer) = Int(x)
_bracket_leaf(x) = throw(ArgumentError("bracket leaves must be integer final-state labels, got $x"))

function _collect_final_labels!(labels::Vector{Int}, tree)
    if tree isa Integer
        push!(labels, Int(tree))
    elseif tree isa Tuple && length(tree) >= 2
        for child in tree
            _collect_final_labels!(labels, child)
        end
    else
        throw(ArgumentError("topology bracket must be a nested tuple of at least two integer leaves"))
    end
    return labels
end

function _assign_lines!(line_ind_for_address::Dict{Any, Int}, next_internal::Base.RefValue{Int}, address)
    if address isa Integer
        line_ind_for_address[address] = _bracket_leaf(address)
        return line_ind_for_address[address]
    end
    address isa Tuple && length(address) >= 2 ||
        throw(ArgumentError("topology bracket must be a nested tuple of at least two integer leaves"))
    for child in address
        _assign_lines!(line_ind_for_address, next_internal, child)
    end
    line = next_internal[]
    next_internal[] += 1
    line_ind_for_address[address] = line
    return line
end

function _collect_vertex_addresses_preorder!(addresses, address)
    address isa Tuple && length(address) >= 2 || return addresses
    push!(addresses, address)
    for child in address
        _collect_vertex_addresses_preorder!(addresses, child)
    end
    return addresses
end

function _count_vertices(tree)
    tree isa Integer && return 0
    tree isa Tuple && length(tree) >= 2 ||
        throw(ArgumentError("topology bracket must be a nested tuple of at least two integer leaves"))
    return 1 + sum(_count_vertices, tree)
end

function _address_final_labels(address)
    labels = _collect_final_labels!(Int[], address)
    return Tuple(labels)
end

function DecayTopology(tree::Tuple)
    final_labels = sort!(_collect_final_labels!(Int[], tree))
    final_labels == collect(1:length(final_labels)) ||
        throw(ArgumentError("bracket leaves must be the consecutive final labels 1:n"))
    nfinal = length(final_labels)
    nvertices = _count_vertices(tree)
    nlines = nfinal + nvertices
    line_ind_for_address = Dict{Any, Int}()
    _assign_lines!(line_ind_for_address, Ref(nfinal + 1), tree)
    vertex_addresses = _collect_vertex_addresses_preorder!(Any[], tree)
    relation = zeros(Int, nlines, nvertices)
    child_order = Vector{Tuple{Vararg{Int}}}(undef, nvertices)
    for (vertex_ind, address) in pairs(vertex_addresses)
        parent = line_ind_for_address[address]
        relation[parent, vertex_ind] = -1
        children = Tuple(line_ind_for_address[child] for child in address)
        for child in children
            relation[child, vertex_ind] = 1
        end
        child_order[vertex_ind] = children
    end
    return _decay_topology_from_relation(relation; root = nlines, finals = Tuple(1:nfinal), child_order = Tuple(child_order))
end

relation(topology::DecayTopology) = topology.relation
root_line_ind(topology::DecayTopology) = topology.root

"""
    final_line_inds(topology)

Return the final-state line ids in canonical order (`1:nfinal`).
"""
final_line_inds(topology::DecayTopology) = topology.finals

nlines(::DecayTopology{Nl}) where {Nl} = Nl
nvertices(::DecayTopology{Nl, Nv}) where {Nl, Nv} = Nv
nfinal(::DecayTopology{Nl, Nv, Nf}) where {Nl, Nv, Nf} = Nf

_line_range(topology::DecayTopology) = Base.OneTo(nlines(topology))
_vertex_range(topology::DecayTopology) = Base.OneTo(nvertices(topology))
_allunique(xs) = length(unique(xs)) == length(xs)

"""
    incoming_line_inds(topology, vertex_ind)

Return incoming graph lines at topology vertex `vertex_ind` in canonical line
order. Valid rooted-tree topologies have exactly one incoming line per vertex;
use [`incoming_line_ind`](@ref) when that invariant is required.
"""
function incoming_line_inds(topology::DecayTopology, vertex_ind::Integer)
    _require_vertex(topology, vertex_ind)
    return [line_ind for line_ind in _line_range(topology) if relation(topology)[line_ind, vertex_ind] == -1]
end

"""
    outgoing_line_inds(topology, vertex_ind)

Return outgoing graph lines at topology vertex `vertex_ind` in canonical line
order as read from [`relation`](@ref). Use [`child_line_inds`](@ref) when the
stored child order matters for traversal or helicity conventions.
"""
function outgoing_line_inds(topology::DecayTopology, vertex_ind::Integer)
    _require_vertex(topology, vertex_ind)
    return [line_ind for line_ind in _line_range(topology) if relation(topology)[line_ind, vertex_ind] == 1]
end

"""
    incoming_line_ind(topology, vertex_ind)

Return the unique incoming graph line at topology vertex `vertex_ind`.
"""
function incoming_line_ind(topology::DecayTopology, vertex_ind::Integer)
    lines = incoming_line_inds(topology, vertex_ind)
    length(lines) == 1 || throw(ArgumentError("vertex_ind $vertex_ind does not have exactly one incoming line"))
    return only(lines)
end

function produced_by(topology::DecayTopology, line_ind::Integer)
    _require_line_ind(topology, line_ind)
    vertex_inds = [
        vertex_ind for vertex_ind in _vertex_range(topology) if relation(topology)[line_ind, vertex_ind] == 1
    ]
    return isempty(vertex_inds) ? nothing : only(vertex_inds)
end

function consumed_by(topology::DecayTopology, line_ind::Integer)
    _require_line_ind(topology, line_ind)
    vertex_inds = [
        vertex_ind for vertex_ind in _vertex_range(topology) if relation(topology)[line_ind, vertex_ind] == -1
    ]
    return isempty(vertex_inds) ? nothing : only(vertex_inds)
end

isroot_line_ind(topology::DecayTopology, line_ind::Integer) = line_ind == root_line_ind(topology)
isfinal_line_ind(topology::DecayTopology, line_ind::Integer) = line_ind in final_line_inds(topology)
isinternal_line_ind(topology::DecayTopology, line_ind::Integer) =
    !isroot_line_ind(topology, line_ind) && !isfinal_line_ind(topology, line_ind)

internal_line_inds(topology::DecayTopology) =
    [line_ind for line_ind in _line_range(topology) if isinternal_line_ind(topology, line_ind)]

has_canonical_line_order(topology::DecayTopology) =
    Tuple(final_line_inds(topology)) == Tuple(1:nfinal(topology)) &&
    internal_line_inds(topology) == collect((nfinal(topology) + 1):(nlines(topology) - 1)) &&
    root_line_ind(topology) == nlines(topology)

function _require_line_ind(topology::DecayTopology, line_ind::Integer)
    line_ind in _line_range(topology) || throw(ArgumentError("line_ind $line_ind is outside 1:$(nlines(topology))"))
    return nothing
end

function _require_vertex(topology::DecayTopology, vertex_ind::Integer)
    vertex_ind in _vertex_range(topology) ||
        throw(ArgumentError("vertex_ind $vertex_ind is outside 1:$(nvertices(topology))"))
    return nothing
end

function _visit_lines!(seen::Set{Int}, stack::Set{Int}, topology::DecayTopology, line_ind::Int)
    line_ind in stack && throw(ArgumentError("topology contains a directed cycle at line_ind $line_ind"))
    push!(stack, line_ind)
    push!(seen, line_ind)
    vertex_ind = consumed_by(topology, line_ind)
    if vertex_ind !== nothing
        for child in outgoing_line_inds(topology, vertex_ind)
            _visit_lines!(seen, stack, topology, child)
        end
    end
    delete!(stack, line_ind)
    return seen
end

"""
    validate_topology(topology)

Validate that a topology is a connected rooted tree with an explicit root line.
Throws `ArgumentError` on invalid input and returns `true` otherwise.
"""
function validate_topology(topology::DecayTopology)
    all(x -> x in (-1, 0, 1), relation(topology)) ||
        throw(ArgumentError("relation entries must be -1, 0, or +1"))

    _require_line_ind(topology, root_line_ind(topology))
    all(line_ind -> line_ind in _line_range(topology), final_line_inds(topology)) ||
        throw(ArgumentError("all final lines must be valid line ids"))
    _allunique(final_line_inds(topology)) || throw(ArgumentError("final lines must be unique"))
    !isfinal_line_ind(topology, root_line_ind(topology)) ||
        throw(ArgumentError("root line cannot also be a final line"))

    nlines(topology) == nfinal(topology) + nvertices(topology) ||
        throw(ArgumentError("tree with explicit root must satisfy nlines == nfinal + nvertices"))

    for vertex_ind in _vertex_range(topology)
        length(incoming_line_inds(topology, vertex_ind)) == 1 ||
            throw(ArgumentError("vertex_ind $vertex_ind must have exactly one incoming line"))
        length(outgoing_line_inds(topology, vertex_ind)) >= 2 ||
            throw(ArgumentError("vertex_ind $vertex_ind must have at least two outgoing lines"))
        sort(collect(child_line_inds(topology, vertex_ind))) == outgoing_line_inds(topology, vertex_ind) ||
            throw(ArgumentError("child_order for vertex_ind $vertex_ind must match outgoing lines"))
    end

    produced_by(topology, root_line_ind(topology)) === nothing ||
        throw(ArgumentError("root line must not be produced by a vertex"))
    consumed_by(topology, root_line_ind(topology)) !== nothing ||
        throw(ArgumentError("root line must be consumed by one vertex"))

    for line_ind in final_line_inds(topology)
        produced_by(topology, line_ind) !== nothing ||
            throw(ArgumentError("final line_ind $line_ind must be produced by one vertex"))
        consumed_by(topology, line_ind) === nothing ||
            throw(ArgumentError("final line_ind $line_ind must not be consumed by a vertex"))
    end

    for line_ind in internal_line_inds(topology)
        produced_by(topology, line_ind) !== nothing ||
            throw(ArgumentError("internal line_ind $line_ind must be produced by one vertex"))
        consumed_by(topology, line_ind) !== nothing ||
            throw(ArgumentError("internal line_ind $line_ind must be consumed by one vertex"))
    end

    seen = _visit_lines!(Set{Int}(), Set{Int}(), topology, root_line_ind(topology))
    length(seen) == nlines(topology) ||
        throw(ArgumentError("topology must be connected to the root line"))

    return true
end

"""
    line_ind_for(topology, address)

Resolve a bracket address such as `(1, 2)` or `((1, 2), 3)` to the internal
line id used by the flat topology. The root address resolves to `root_line_ind`.
"""
function line_ind_for(topology::DecayTopology, address)
    labels = _address_final_labels(address)
    for line_ind in _line_range(topology)
        _final_descendants_tuple(topology, line_ind) == labels && return Int(line_ind)
    end
    throw(ArgumentError("address $address does not correspond to a line_ind in this topology"))
end

"""
    vertex_ind_for(topology, address)

Resolve the parent bracket address of a binary decay to the vertex id. For
example, `(1, 2)` resolves to the vertex `(1,2) -> 1,2`.
"""
function vertex_ind_for(topology::DecayTopology, address)
    line_ind = line_ind_for(topology, address)
    vertex_ind = consumed_by(topology, line_ind)
    vertex_ind === nothing && throw(ArgumentError("address $address does not correspond to a decay vertex"))
    return vertex_ind
end

"""
    vertex_address(topology, vertex_ind)

Return the bracket address of a decay vertex, i.e. the same key used in
[`vertex_ind_for`](@ref) and in [`DecayChain`](@ref) vertex payloads.
"""
function vertex_address(topology::DecayTopology, vertex_ind::Integer)
    _require_vertex(topology, vertex_ind)
    return _line_ind_address(topology, incoming_line_ind(topology, vertex_ind))
end

function _compact_vertex_label(address)
    address isa Integer && return string(address)
    return "(" * join((_compact_vertex_label(child) for child in address), ",") * ")"
end

function _line_ind_address(topology::DecayTopology, line_ind::Integer)
    _require_line_ind(topology, line_ind)
    isfinal_line_ind(topology, line_ind) && return Int(line_ind)
    children = child_line_inds(topology, consumed_by(topology, line_ind))
    return Tuple(_line_ind_address(topology, child) for child in children)
end

function final_descendants(topology::DecayTopology, line_ind::Integer)
    _require_line_ind(topology, line_ind)
    vertex_ind = consumed_by(topology, line_ind)
    vertex_ind === nothing && return [Int(line_ind)]
    result = Int[]
    for child in child_line_inds(topology, vertex_ind)
        append!(result, final_descendants(topology, child))
    end
    return result
end

"""
    child_line_inds(topology, vertex_ind)

Return the outgoing lines of topology vertex `vertex_ind` in stored child order.
"""
function child_line_inds(topology::DecayTopology, vertex_ind::Integer)
    _require_vertex(topology, vertex_ind)
    children = topology.child_order[Int(vertex_ind)]
    return SVector{length(children), Int}(children)
end

"""
    arity(topology, vertex_ind)

Return the number of ordered child lines attached to topology vertex
`vertex_ind`.
"""
arity(topology::DecayTopology, vertex_ind::Integer) = length(child_line_inds(topology, vertex_ind))

"""
    is_binary_vertex(topology, vertex_ind)

Return `true` when topology vertex `vertex_ind` has exactly two outgoing child
lines.
"""
is_binary_vertex(topology::DecayTopology, vertex_ind::Integer) = arity(topology, vertex_ind) == 2

"""
    line_inds_at_vertex(topology, vertex_ind)

Return `(parent, children...)` for a topology vertex. Unlike
[`vertex_line_inds`](@ref), this helper supports vertices with any arity.
"""
function line_inds_at_vertex(topology::DecayTopology, vertex_ind::Integer)
    parent = incoming_line_ind(topology, vertex_ind)
    children = child_line_inds(topology, vertex_ind)
    return (parent, Tuple(children)...)
end

"""
    vertex_line_inds(topology, vertex_ind)

Return `(parent, child1, child2)` for a binary vertex. The children are returned
in stored child order.
"""
function vertex_line_inds(topology::DecayTopology, vertex_ind::Integer)
    parent = incoming_line_ind(topology, vertex_ind)
    children = child_line_inds(topology, vertex_ind)
    length(children) == 2 ||
        throw(ArgumentError("vertex_ind $vertex_ind is not binary"))
    return (parent, children[1], children[2])
end

function _line_ind_label(line_ind::Int, labels)
    labels === nothing && return string(line_ind)
    return string(labels[line_ind])
end

function _bracket_for_line_ind(topology::DecayTopology, line_ind::Int, labels)
    vertex_ind = consumed_by(topology, line_ind)
    vertex_ind === nothing && return _line_ind_label(line_ind, labels)
    pieces = [_bracket_for_line_ind(topology, child, labels) for child in child_line_inds(topology, vertex_ind)]
    return "(" * join(pieces, ",") * ")"
end

"""
    bracket(topology; labels=nothing)

Reconstruct a compact bracket notation from the flat topology.
"""
bracket(topology::DecayTopology; labels = nothing) =
    _bracket_for_line_ind(topology, root_line_ind(topology), labels)

"""
    bracket_notation(topology; labels=nothing)

Return the compact nested-tuple notation for a topology, for example
`"((1,2),3)"`.
"""
bracket_notation(topology::DecayTopology; labels = nothing) =
    bracket(topology; labels)

function _final_descendants_tuple(topology::DecayTopology, line_ind::Integer)
    _require_line_ind(topology, line_ind)
    vertex_ind = consumed_by(topology, line_ind)
    vertex_ind === nothing && return (Int(line_ind),)
    labels = Int[]
    for child in child_line_inds(topology, vertex_ind)
        append!(labels, _final_descendants_tuple(topology, child))
    end
    return Tuple(labels)
end

_indices_for_line_ind(topology::DecayTopology, line_ind::Integer) =
    _final_descendants_tuple(topology, line_ind)

function _subtree_line_inds(topology::DecayTopology, line_ind::Integer)
    line_inds = Int[Int(line_ind)]
    vertex_ind = consumed_by(topology, line_ind)
    vertex_ind === nothing && return line_inds
    for child in child_line_inds(topology, vertex_ind)
        append!(line_inds, _subtree_line_inds(topology, child))
    end
    return line_inds
end

function _child_containing_line_ind(topology::DecayTopology, vertex_ind::Integer, line_ind::Integer)
    for child in child_line_inds(topology, vertex_ind)
        if line_ind in _subtree_line_inds(topology, child)
            return child
        end
    end
    throw(ArgumentError("line_ind $line_ind is not below vertex $vertex_ind"))
end

function _root_vertex_ind(topology::DecayTopology)
    vertex_ind = consumed_by(topology, root_line_ind(topology))
    vertex_ind === nothing && throw(ArgumentError("topology root line is not consumed by any vertex"))
    return vertex_ind
end
