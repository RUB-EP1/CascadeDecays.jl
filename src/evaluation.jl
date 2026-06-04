_indices_for_line(topology::DecayTopology, line::Integer) =
    Tuple(final_descendants(topology, line))

_neg_indices(indices::Tuple) = Tuple(-i for i in indices)

function _child_containing_line(topology::DecayTopology, vertex::Integer, line::Integer)
    for child in child_lines(topology, vertex)
        descendants = final_descendants(topology, child)
        if child == line || line in descendants || !isfinalline(topology, line) && !isrootline(topology, line) && line in _subtree_lines(topology, child)
            return child
        end
    end
    throw(ArgumentError("line $line is not below vertex $vertex"))
end

function _subtree_lines(topology::DecayTopology, line::Integer)
    lines = Int[Int(line)]
    vertex = consumed_by(topology, line)
    vertex === nothing && return lines
    for child in child_lines(topology, vertex)
        append!(lines, _subtree_lines(topology, child))
    end
    return lines
end

function _root_vertex(topology::DecayTopology)
    vertex = consumed_by(topology, rootline(topology))
    vertex === nothing && throw(ArgumentError("topology root line is not consumed by any vertex"))
    return vertex
end

abstract type AbstractInitialFrame end

"""
    HelicityRootFrame()

Start generated helicity-angle programs by transforming to the root/system
helicity frame. This is the default for fully general four-vectors.
"""
struct HelicityRootFrame <: AbstractInitialFrame end

"""
    CurrentFrame()

Start generated helicity-angle programs in the current axes. Use this when the
input four-vectors are already expressed in the intended system frame.
"""
struct CurrentFrame <: AbstractInitialFrame end

_initial_frame_program(topology::DecayTopology, ::HelicityRootFrame) =
    (ToHelicityFrame(_indices_for_line(topology, rootline(topology))),)
_initial_frame_program(::DecayTopology, ::CurrentFrame) = ()

"""
    helicity_angle_program(topology, vertex; initial_frame=HelicityRootFrame())

Build an `InstructionalDecayTrees.jl` instruction program that measures the
local `(cosθ, ϕ)` angle for `vertex` in helicity convention.
"""
function helicity_angle_program(
    topology::DecayTopology,
    target_vertex::Integer;
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
)
    return helicity_angle_program(topology, Val(Int(target_vertex)); initial_frame)
end

function helicity_angle_program(
    topology::DecayTopology,
    ::Val{target_vertex};
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
) where {target_vertex}
    target_vertex in Base.OneTo(nvertices(topology)) ||
        throw(ArgumentError("vertex $target_vertex is outside 1:$(nvertices(topology))"))

    target_parent = incoming_line(topology, target_vertex)
    program = _initial_frame_program(topology, initial_frame)

    current_vertex = _root_vertex(topology)

    while true
        parent = incoming_line(topology, current_vertex)
        children = child_lines(topology, current_vertex)
        if parent == target_parent
            return (
                program...,
                MeasureCosThetaPhi(Symbol(:v, target_vertex), _indices_for_line(topology, children[1])),
            )
        end
        next_child = _child_containing_line(topology, current_vertex, target_parent)
        program = (program..., ToHelicityFrame(_indices_for_line(topology, next_child)))
        current_vertex = consumed_by(topology, next_child)
        current_vertex === nothing &&
            throw(ArgumentError("target vertex $target_vertex is not reachable from the root"))
    end
end

"""
    helicity_angle_programs(topology; initial_frame=HelicityRootFrame())

Return one helicity-angle measurement program per topology vertex.
"""
helicity_angle_programs(
    topology::DecayTopology;
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
) =
    ntuple(v -> helicity_angle_program(topology, Val(v); initial_frame), nvertices(topology))

function _sum_objects(objs, indices::Tuple)
    return reduce(+, (objs[i] for i in indices))
end

"""
    cascade_kinematics(topology, system, objs; initial_frame=HelicityRootFrame())

Compute a `CascadeKinematics` input from final-state four-vectors `objs`.
Internal invariant masses are derived from final-state descendants, and vertex
angles are computed with generated helicity-angle programs.
"""
function cascade_kinematics(
    topology::DecayTopology,
    system::CascadeSystem,
    objs;
    initial_frame::AbstractInitialFrame=HelicityRootFrame(),
)
    internal_masses2 = Tuple(mass(_sum_objects(objs, _indices_for_line(topology, line)))^2 for line in internal_lines(topology))
    programs = helicity_angle_programs(topology; initial_frame)
    angle_results = map(programs) do program
        _, result = apply_decay_instruction(program, objs)
        only(values(result))
    end
    return CascadeKinematics(topology, system; internal_masses2, vertex_angles=Tuple(angle_results))
end

function routed_vertex_amplitude(vertex, masses2, helicities, spins, angles)
    throw(MethodError(routed_vertex_amplitude, (vertex, masses2, helicities, spins, angles)))
end

function routed_vertex_amplitude(
    vertex::ThreeBodyDecays.VertexFunction,
    masses2,
    helicities,
    spins,
    angles,
)
    two_j0, _, _ = spins
    two_λ0, two_λ1, two_λ2 = helicities
    recoupling =
        _particle_two_phase(spins[3], two_λ2) *
        ThreeBodyDecays.amplitude(vertex.h, (two_λ1, two_λ2), spins)
    formfactor = vertex.ff(masses2...)
    two_Δλ = two_λ1 - two_λ2
    rotation = conj(ThreeBodyDecays.wignerD_doublearg(two_j0, two_λ0, two_Δλ, angles.ϕ, angles.cosθ, 0))
    return rotation * recoupling * formfactor
end

function _particle_two_phase(two_j2::Integer, two_λ2::Integer)
    exponent_num = two_j2 - two_λ2
    iseven(exponent_num) ||
        throw(ArgumentError("particle-2 phase requires two_j2 - two_λ2 to be even"))
    return isodd(div(exponent_num, 2)) ? -1 : 1
end

function routed_vertex_amplitude(chain::DecayChain, system::CascadeSystem, x::CascadeKinematics, two_λs, vertex::Integer)
    masses2 = vertex_masses2(chain, x, vertex)
    helicities = vertex_helicities(chain, two_λs, vertex)
    spins = vertex_spins(chain, system, vertex)
    angles = vertex_angles(x, vertex)
    return routed_vertex_amplitude(chain.vertices[vertex], masses2, helicities, spins, angles)
end

function routed_propagator_product(chain::DecayChain, x::CascadeKinematics)
    return prod(zip(chain.propagators, propagating_lines(chain))) do (propagator, line)
        propagator(line_invariant(x, line))
    end
end

# helicity-axis helpers (same indexing as ThreeBodyDecays: div(two_j + two_λ, 2) + 1)
_helicity_axis_length(two_j::Integer) = two_j + 1
_helicity_index(two_λ::Integer, two_j::Integer) = div(two_j + two_λ, 2) + 1

function _external_amplitude_indices(
    chain::DecayChain,
    system::CascadeSystem,
    external_two_λs::SystemSpins,
)
    two_js = line_two_js(chain, system)
    final_indices = ntuple(
        i -> _helicity_index(external_two_λs.finals[i], two_js[finallines(chain)[i]]),
        nfinal(chain),
    )
    root_index = _helicity_index(external_two_λs.two_h0, root_two_j(system))
    return (final_indices..., root_index)
end

"""
    _vertex_factor(chain, system, x, vertex)

Local vertex amplitude ``V_{λ_0 λ_1 λ_2}`` on the three lines of `vertex`,
as a dense array (cf. `VRk` / `Vij` in `ThreeBodyDecays.aligned_amplitude`).
"""
function _vertex_factor(
    chain::DecayChain,
    system::CascadeSystem,
    x::CascadeKinematics,
    vertex::Integer,
)
    l0, l1, l2 = vertex_lines(chain, vertex)
    two_j0, two_j1, two_j2 = vertex_spins(chain, system, vertex)
    masses2 = vertex_masses2(chain, x, vertex)
    spins = (two_j0, two_j1, two_j2)
    angles = vertex_angles(x, vertex)
    vertex_payload = chain.vertices[vertex]
    V = [
        routed_vertex_amplitude(vertex_payload, masses2, (two_λ0, two_λ1, two_λ2), spins, angles)
        for two_λ0 in (-two_j0):2:two_j0, two_λ1 in (-two_j1):2:two_j1, two_λ2 in (-two_j2):2:two_j2
    ]
    return reshape(
        V,
        _helicity_axis_length(two_j0),
        _helicity_axis_length(two_j1),
        _helicity_axis_length(two_j2),
    ),
    (l0, l1, l2)
end

function _multiply_vertex_into_lines!(
    F::AbstractArray{T,N},
    V::AbstractArray,
    lines::NTuple{3,Int},
) where {T,N}
    l0, l1, l2 = lines
    expand_sizes = ntuple(
        d -> d == l0 ? size(V, 1) : d == l1 ? size(V, 2) : d == l2 ? size(V, 3) : 1,
        Val(N),
    )
    F .*= reshape(V, expand_sizes)
    return F
end

"""Element type for line / external amplitude arrays (cf. `zeros(typeof(lineshape), …)` in ThreeBodyDecays)."""
function _line_amplitude_eltype(chain::DecayChain, system::CascadeSystem, x::CascadeKinematics)
    P = routed_propagator_product(chain, x)
    v = 1
    two_j0, two_j1, two_j2 = vertex_spins(chain, system, v)
    helicities = (
        first((-two_j0):2:two_j0),
        first((-two_j1):2:two_j1),
        first((-two_j2):2:two_j2),
    )
    V0 = routed_vertex_amplitude(
        chain.vertices[v],
        vertex_masses2(chain, x, v),
        helicities,
        (two_j0, two_j1, two_j2),
        vertex_angles(x, v),
    )
    return promote_type(typeof(V0), typeof(P))
end

"""
    line_amplitude_tensor(chain, system, x)

Line-indexed amplitude buffer ``F_{λ_{\\mathrm{line}_1} \\ldots} =
\\prod_v V_v`` before summing internal propagator helicities.
"""
function line_amplitude_tensor(
    chain::DecayChain,
    system::CascadeSystem,
    x::CascadeKinematics,
)
    T = _line_amplitude_eltype(chain, system, x)
    two_js = line_two_js(chain, system)
    line_sizes = ntuple(line -> _helicity_axis_length(two_js[line]), nlines(chain))
    F = ones(T, line_sizes...)
    for vertex in 1:nvertices(chain)
        V, lines = _vertex_factor(chain, system, x, vertex)
        _multiply_vertex_into_lines!(F, V, lines)
    end
    return F
end

function _permute_external(F::AbstractArray{T,N}, ext_dims::NTuple{Ne,Int}) where {T,N,Ne}
    Ne == N && return F
    other_dims = Tuple(d for d in 1:N if d ∉ ext_dims)
    Fp = permutedims(F, (ext_dims..., other_dims...))
    return dropdims(Fp; dims=ntuple(i -> Ne + i, Val(N - Ne)))
end

# Static `@tullio` methods (generated at load time), cf. `ThreeBodyDecays.amplitude(dc, σs)`.
for Ne in 1:8, Ni in 1:4
    ext_syms = Tuple(Symbol(:_λ, k) for k in 1:Ne)
    int_syms = Tuple(Symbol(:_m, k) for k in 1:Ni)
    fp_syms = (ext_syms..., int_syms...)
    A_shape = Tuple(:(size(Fp, $k)) for k in 1:Ne)
    N = Ne + Ni
    @eval function _tullio_sum_internal(Fp::AbstractArray{Ta,$N}, ::Val{$Ne}, ::Val{$Ni}) where {Ta}
        A = zeros(Ta, $(A_shape...))
        @tullio $(Expr(:ref, :A, ext_syms...)) := $(Expr(:ref, :Fp, fp_syms...))
        return A
    end
end

function _tullio_sum_internal(
    Fp::AbstractArray{Ta,N},
    ::Val{Ne},
    ::Val{Ni},
) where {Ta,N,Ne,Ni}
    int_dims = ntuple(i -> Ne + i, Val(Ni))
    return dropdims(sum(Fp; dims=int_dims), dims=int_dims)
end

"""
    external_helicity_amplitude(F, ext_dims, int_dims)

Sum internal propagator helicities and return the external-helicity array.
`ext_dims` / `int_dims` are line ids (finals + root, then propagating lines).
"""
function external_helicity_amplitude(
    F::AbstractArray,
    ext_dims::Tuple{Int,Vararg{Int}},
    int_dims::Tuple{Int,Vararg{Int}},
)
    isempty(int_dims) && return _permute_external(F, ext_dims)
    Fp = permutedims(F, (ext_dims..., int_dims...))
    return _tullio_sum_internal(Fp, Val(length(ext_dims)), Val(length(int_dims)))
end

"""
    amplitude(chain, system, x)

Route masses, angles, and spins through the cascade and return the full amplitude
array in the external helicity space. Final-state axes follow [`finallines`](@ref)
order; the root helicity is the last axis.

The evaluation follows `ThreeBodyDecays.jl`: build vertex factors, multiply into
a line-indexed buffer, then contract internal helicities with `@tullio`.
"""
function amplitude(
    chain::DecayChain{Nf,Np},
    system::CascadeSystem,
    x::CascadeKinematics,
) where {Nf,Np}
    P_prod = routed_propagator_product(chain, x)
    F = line_amplitude_tensor(chain, system, x)
    ext_dims = (Tuple(finallines(chain))..., rootline(chain))
    int_dims = Tuple(propagating_lines(chain))
    A = external_helicity_amplitude(F, ext_dims, int_dims)
    return P_prod * A
end

"""
    amplitude(chain, system, x, external_two_λs)

Return one helicity component of [`amplitude`](@ref)(`chain`, `system`, `x`).
`external_two_λs` is a [`SystemHelicities`](@ref) value (alias of [`SystemSpins`](@ref):
positional finals, root via `two_h0=...` or `h0=...`).
"""
function amplitude(
    chain::DecayChain,
    system::CascadeSystem,
    x::CascadeKinematics,
    external_two_λs::SystemSpins,
)
    A = amplitude(chain, system, x)
    return A[_external_amplitude_indices(chain, system, external_two_λs)...]
end
