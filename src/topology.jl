"""
    DecayTopology(relation; root, finals)

Flat line-vertex incidence graph for a connected binary cascade tree.

Rows of `relation` are lines and columns are vertices. The sign convention is:

- `-1`: line enters a vertex
- `+1`: line leaves a vertex
- `0`: line is unrelated to the vertex

The root/mother line is included explicitly. Final-state lines are listed in
`finals`.
"""
struct DecayTopology{Nl,Nv,Nf,T<:Integer,L}
    relation::SMatrix{Nl,Nv,T,L}
    root::Int
    finals::SVector{Nf,Int}
end

function DecayTopology(
    relation::SMatrix{Nl,Nv,T,L},
    root::Integer,
    finals::SVector{Nf,Int};
    validate::Bool = true,
) where {Nl,Nv,Nf,T<:Integer,L}
    topology = DecayTopology{Nl,Nv,Nf,T,L}(relation, Int(root), finals)
    validate && validate_topology(topology)
    return topology
end

function DecayTopology(
    relation::AbstractMatrix{T};
    root::Integer,
    finals,
    validate::Bool = true,
) where {T<:Integer}
    Nl, Nv = size(relation)
    final_tuple = Tuple(Int(line) for line in finals)
    final_lines = SVector{length(final_tuple),Int}(final_tuple)
    static_relation = SMatrix{Nl,Nv,T,Nl * Nv}(relation)
    return DecayTopology(static_relation, root, final_lines; validate)
end

_bracket_leaf(x::Integer) = Int(x)
_bracket_leaf(x) = throw(ArgumentError("bracket leaves must be integer final-state labels, got $x"))

function _collect_final_labels!(labels::Vector{Int}, tree)
    if tree isa Integer
        push!(labels, Int(tree))
    elseif tree isa Tuple && length(tree) == 2
        _collect_final_labels!(labels, tree[1])
        _collect_final_labels!(labels, tree[2])
    else
        throw(ArgumentError("topology bracket must be a binary nested tuple of integer leaves"))
    end
    return labels
end

function _assign_lines!(line_for_node::Dict{Any,Int}, next_internal::Base.RefValue{Int}, tree)
    if tree isa Integer
        line_for_node[tree] = _bracket_leaf(tree)
        return line_for_node[tree]
    end
    tree isa Tuple && length(tree) == 2 ||
        throw(ArgumentError("topology bracket must be a binary nested tuple of integer leaves"))
    _assign_lines!(line_for_node, next_internal, tree[1])
    _assign_lines!(line_for_node, next_internal, tree[2])
    line = next_internal[]
    next_internal[] += 1
    line_for_node[tree] = line
    return line
end

function _collect_vertices_preorder!(vertices, tree)
    tree isa Tuple && length(tree) == 2 || return vertices
    push!(vertices, tree)
    _collect_vertices_preorder!(vertices, tree[1])
    _collect_vertices_preorder!(vertices, tree[2])
    return vertices
end

function _address_final_labels(address)
    labels = _collect_final_labels!(Int[], address)
    return Tuple(labels)
end

"""
    DecayTopology(bracket)

Construct a flat line-vertex incidence topology from binary bracket notation,
for example `DecayTopology((((1, 2), 3), 4))`.

The generated convention is final-state lines first, internal lines next in
postorder, and the root/mother line last. Vertex columns are root-first.
"""
function DecayTopology(tree::Tuple)
    final_labels = sort!(_collect_final_labels!(Int[], tree))
    final_labels == collect(1:length(final_labels)) ||
        throw(ArgumentError("bracket leaves must be the consecutive final labels 1:n"))
    nfinal = length(final_labels)
    nvertices = nfinal - 1
    nlines = nfinal + nvertices
    line_for_node = Dict{Any,Int}()
    _assign_lines!(line_for_node, Ref(nfinal + 1), tree)
    vertices = _collect_vertices_preorder!(Any[], tree)
    relation = zeros(Int, nlines, nvertices)
    for (vertex_ind, node) in pairs(vertices)
        parent = line_for_node[node]
        child1 = line_for_node[node[1]]
        child2 = line_for_node[node[2]]
        relation[parent, vertex_ind] = -1
        relation[child1, vertex_ind] = 1
        relation[child2, vertex_ind] = 1
    end
    return DecayTopology(relation; root = nlines, finals = Tuple(1:nfinal))
end

relation(topology::DecayTopology) = topology.relation
rootline(topology::DecayTopology) = topology.root
finallines(topology::DecayTopology) = topology.finals
nlines(::DecayTopology{Nl}) where {Nl} = Nl
nvertices(::DecayTopology{Nl,Nv}) where {Nl,Nv} = Nv
nfinal(::DecayTopology{Nl,Nv,Nf}) where {Nl,Nv,Nf} = Nf

_line_range(topology::DecayTopology) = Base.OneTo(nlines(topology))
_vertex_range(topology::DecayTopology) = Base.OneTo(nvertices(topology))
_allunique(xs) = length(unique(xs)) == length(xs)

function incoming_lines(topology::DecayTopology, vertex_ind::Integer)
    _require_vertex(topology, vertex_ind)
    return [line for line in _line_range(topology) if relation(topology)[line, vertex_ind] == -1]
end

function outgoing_lines(topology::DecayTopology, vertex_ind::Integer)
    _require_vertex(topology, vertex_ind)
    return [line for line in _line_range(topology) if relation(topology)[line, vertex_ind] == 1]
end

function incoming_line(topology::DecayTopology, vertex_ind::Integer)
    lines = incoming_lines(topology, vertex_ind)
    length(lines) == 1 || throw(ArgumentError("vertex_ind $vertex_ind does not have exactly one incoming line"))
    return only(lines)
end

function produced_by(topology::DecayTopology, line::Integer)
    _require_line(topology, line)
    vertex_inds = [
        vertex_ind for vertex_ind in _vertex_range(topology) if relation(topology)[line, vertex_ind] == 1
    ]
    return isempty(vertex_inds) ? nothing : only(vertex_inds)
end

function consumed_by(topology::DecayTopology, line::Integer)
    _require_line(topology, line)
    vertex_inds = [
        vertex_ind for vertex_ind in _vertex_range(topology) if relation(topology)[line, vertex_ind] == -1
    ]
    return isempty(vertex_inds) ? nothing : only(vertex_inds)
end

isrootline(topology::DecayTopology, line::Integer) = line == rootline(topology)
isfinalline(topology::DecayTopology, line::Integer) = line in finallines(topology)
isinternalline(topology::DecayTopology, line::Integer) =
    !isrootline(topology, line) && !isfinalline(topology, line)

internal_lines(topology::DecayTopology) =
    [line for line in _line_range(topology) if isinternalline(topology, line)]

has_canonical_line_order(topology::DecayTopology) =
    Tuple(finallines(topology)) == Tuple(1:nfinal(topology)) &&
    internal_lines(topology) == collect((nfinal(topology) + 1):(nlines(topology) - 1)) &&
    rootline(topology) == nlines(topology)

function _require_line(topology::DecayTopology, line::Integer)
    line in _line_range(topology) || throw(ArgumentError("line $line is outside 1:$(nlines(topology))"))
    return nothing
end

function _require_vertex(topology::DecayTopology, vertex_ind::Integer)
    vertex_ind in _vertex_range(topology) ||
        throw(ArgumentError("vertex_ind $vertex_ind is outside 1:$(nvertices(topology))"))
    return nothing
end

function _visit_lines!(seen::Set{Int}, stack::Set{Int}, topology::DecayTopology, line::Int)
    line in stack && throw(ArgumentError("topology contains a directed cycle at line $line"))
    push!(stack, line)
    push!(seen, line)
    vertex_ind = consumed_by(topology, line)
    if vertex_ind !== nothing
        for child in outgoing_lines(topology, vertex_ind)
            _visit_lines!(seen, stack, topology, child)
        end
    end
    delete!(stack, line)
    return seen
end

"""
    validate_topology(topology)

Validate that a topology is a connected binary tree with an explicit root line.
Throws `ArgumentError` on invalid input and returns `true` otherwise.
"""
function validate_topology(topology::DecayTopology)
    all(x -> x in (-1, 0, 1), relation(topology)) ||
        throw(ArgumentError("relation entries must be -1, 0, or +1"))

    _require_line(topology, rootline(topology))
    all(line -> line in _line_range(topology), finallines(topology)) ||
        throw(ArgumentError("all final lines must be valid line ids"))
    _allunique(finallines(topology)) || throw(ArgumentError("final lines must be unique"))
    !isfinalline(topology, rootline(topology)) ||
        throw(ArgumentError("root line cannot also be a final line"))

    nlines(topology) == nfinal(topology) + nvertices(topology) ||
        throw(ArgumentError("binary tree with explicit root must satisfy nlines == nfinal + nvertices"))

    for vertex_ind in _vertex_range(topology)
        length(incoming_lines(topology, vertex_ind)) == 1 ||
            throw(ArgumentError("vertex_ind $vertex_ind must have exactly one incoming line"))
        length(outgoing_lines(topology, vertex_ind)) == 2 ||
            throw(ArgumentError("vertex_ind $vertex_ind must have exactly two outgoing lines"))
    end

    produced_by(topology, rootline(topology)) === nothing ||
        throw(ArgumentError("root line must not be produced by a vertex"))
    consumed_by(topology, rootline(topology)) !== nothing ||
        throw(ArgumentError("root line must be consumed by one vertex"))

    for line in finallines(topology)
        produced_by(topology, line) !== nothing ||
            throw(ArgumentError("final line $line must be produced by one vertex"))
        consumed_by(topology, line) === nothing ||
            throw(ArgumentError("final line $line must not be consumed by a vertex"))
    end

    for line in internal_lines(topology)
        produced_by(topology, line) !== nothing ||
            throw(ArgumentError("internal line $line must be produced by one vertex"))
        consumed_by(topology, line) !== nothing ||
            throw(ArgumentError("internal line $line must be consumed by one vertex"))
    end

    seen = _visit_lines!(Set{Int}(), Set{Int}(), topology, rootline(topology))
    length(seen) == nlines(topology) ||
        throw(ArgumentError("topology must be connected to the root line"))

    return true
end

function _final_descendants_unordered(topology::DecayTopology, line::Integer)
    _require_line(topology, line)
    vertex_ind = consumed_by(topology, line)
    vertex_ind === nothing && return [Int(line)]
    result = Int[]
    for child in outgoing_lines(topology, vertex_ind)
        append!(result, _final_descendants_unordered(topology, child))
    end
    return result
end

"""
    line_for(topology, address)

Resolve a bracket address such as `(1, 2)` or `((1, 2), 3)` to the internal
line id used by the flat topology. The root address resolves to `rootline`.
"""
function line_for(topology::DecayTopology, address)
    labels = _address_final_labels(address)
    for line in _line_range(topology)
        Tuple(final_descendants(topology, line)) == labels && return Int(line)
    end
    throw(ArgumentError("address $address does not correspond to a line in this topology"))
end

"""
    vertex_for(topology, address)

Resolve the parent bracket address of a binary decay to the vertex id. For
example, `(1, 2)` resolves to the vertex `(1,2) -> 1,2`.
"""
function vertex_for(topology::DecayTopology, address)
    line = line_for(topology, address)
    vertex_ind = consumed_by(topology, line)
    vertex_ind === nothing && throw(ArgumentError("address $address does not correspond to a decay vertex"))
    return vertex_ind
end

"""
    vertex_address(topology, vertex_ind)

Return the bracket address of a decay vertex, i.e. the same key used in
[`vertex_for`](@ref) and in [`DecayChain`](@ref) vertex payloads.
"""
function vertex_address(topology::DecayTopology, vertex_ind::Integer)
    _require_vertex(topology, vertex_ind)
    return _line_address(topology, incoming_line(topology, vertex_ind))
end

function _line_address(topology::DecayTopology, line::Integer)
    _require_line(topology, line)
    isfinalline(topology, line) && return Int(line)
    children = child_lines(topology, consumed_by(topology, line))
    return (_line_address(topology, children[1]), _line_address(topology, children[2]))
end

function _ordered_children(topology::DecayTopology, vertex_ind::Int)
    children = outgoing_lines(topology, vertex_ind)
    return sort(children; by = line -> minimum(_final_descendants_unordered(topology, line)))
end

function final_descendants(topology::DecayTopology, line::Integer)
    _require_line(topology, line)
    vertex_ind = consumed_by(topology, line)
    vertex_ind === nothing && return [Int(line)]
    result = Int[]
    for child in child_lines(topology, vertex_ind)
        append!(result, final_descendants(topology, child))
    end
    return result
end

"""
    child_lines(topology, vertex_ind)

Return the two outgoing lines of topology vertex `vertex_ind` in canonical child
order. This is the ordered pair that should be used for local two-body arguments.
"""
function child_lines(topology::DecayTopology, vertex_ind::Integer)
    _require_vertex(topology, vertex_ind)
    return _ordered_children(topology, Int(vertex_ind))
end

"""
    vertex_lines(topology, vertex_ind)

Return `(parent, child1, child2)` for a binary vertex. The children are ordered
canonically by their final-state descendants.
"""
function vertex_lines(topology::DecayTopology, vertex_ind::Integer)
    parent = incoming_line(topology, vertex_ind)
    children = child_lines(topology, vertex_ind)
    length(children) == 2 ||
        throw(ArgumentError("vertex_ind $vertex_ind does not have exactly two outgoing lines"))
    return (parent, children[1], children[2])
end

function _line_label(line::Int, labels)
    labels === nothing && return string(line)
    return string(labels[line])
end

function _bracket_for_line(topology::DecayTopology, line::Int, labels)
    vertex_ind = consumed_by(topology, line)
    vertex_ind === nothing && return _line_label(line, labels)
    pieces = [_bracket_for_line(topology, child, labels) for child in child_lines(topology, vertex_ind)]
    return "(" * join(pieces, ",") * ")"
end

"""
    bracket(topology; labels=nothing)

Reconstruct a compact bracket notation from the canonical flat topology.
"""
bracket(topology::DecayTopology; labels = nothing) =
    _bracket_for_line(topology, rootline(topology), labels)
