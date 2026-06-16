"""
Three-body resonance–spectator conventions for **tests, benchmarks, and docs**
cross-checking [`ThreeBodyDecays.jl`](https://github.com/mmikhasenko/ThreeBodyDecays.jl).

This module is not loaded by `CascadeDecays` itself. Use it from `test/`,
`benchmark/`, and `docs/` when comparing against `amplitude(dc, po, σs; refζs)`.
"""
module ThreeBodyCompat

using FourVectors
using CascadeDecays
using CascadeDecays:
    DecayTopology,
    apply_decay_instruction,
    child_line_inds,
    consumed_by,
    final_descendants,
    helicity_angle_programs,
    internal_line_inds,
    isfinal_line_ind,
    nfinal,
    root_line_ind,
    vertex_address,
    AbstractInitialFrame,
    CurrentFrame

export three_body_topology,
    three_body_k,
    three_body_spectator,
    three_body_isobar_pair,
    three_body_root_vertex,
    three_body_refζs,
    three_body_σs,
    reference_plane_orientation

"""
    three_body_topology(k)

Return the resonance–spectator [`DecayTopology`](@ref) for a three-body chain
with isobar index `k` in `ThreeBodyDecays.jl`.
"""
function three_body_topology(k::Integer)
    k in (1, 2, 3) || throw(ArgumentError("three_body_topology requires k in (1, 2, 3)"))
    return DecayTopology(_three_body_tree(k))
end

three_body_spectator(k::Integer) = k
three_body_isobar_pair(k::Integer) = k == 1 ? (2, 3) : k == 2 ? (3, 1) : (1, 2)
three_body_root_vertex(k::Integer) = three_body_topology(k) |> t -> vertex_address(t, 1)

function three_body_refζs(k::Integer)
    k in (1, 2, 3) || throw(ArgumentError("three_body_refζs requires k in (1, 2, 3)"))
    return k == 1 ? ntuple(_ -> 1, 4) : ntuple(_ -> k, 4)
end

three_body_refζs(topology::DecayTopology) = three_body_refζs(three_body_k(topology))
_three_body_tree(k::Integer) = (three_body_isobar_pair(k), three_body_spectator(k))

function three_body_k(topology::DecayTopology)
    nfinal(topology) == 3 ||
        throw(ArgumentError("three_body_k requires a three-final topology"))
    root_vertex = consumed_by(topology, root_line_ind(topology))
    for child in child_line_inds(topology, root_vertex)
        if isfinal_line_ind(topology, child)
            return only(final_descendants(topology, child))
        end
    end
    throw(ArgumentError("could not identify three-body spectator line"))
end

function three_body_σs(objs)
    return (
        σ1 = mass(objs[2] + objs[3])^2,
        σ2 = mass(objs[1] + objs[3])^2,
        σ3 = mass(objs[1] + objs[2])^2,
    )
end

_spatial_norm2(p) = p.px^2 + p.py^2 + p.pz^2
_rest_scale2(p) = max(abs(p.E), one(float(abs(p.E))))^2
_effectively_at_rest(p; rtol = 1.0e-12) = _spatial_norm2(p) <= rtol^2 * _rest_scale2(p)

function _sum_line_momentum(topology::DecayTopology, objs, line_ind::Integer)
    return sum(i -> objs[i], final_descendants(topology, line_ind))
end

_root_isobar_line(topology::DecayTopology) = only(internal_line_inds(topology))

function reference_plane_orientation(
        reference_topology::DecayTopology,
        objs;
        initial_frame::AbstractInitialFrame = CurrentFrame(),
    )
    nfinal(reference_topology) == 3 ||
        throw(ArgumentError("reference plane orientation requires a three-final topology"))
    return _reference_plane_orientation(reference_topology, objs, initial_frame)
end

function _reference_plane_orientation(
        reference_topology::DecayTopology,
        objs,
        initial_frame::HelicityRootFrame,
    )
    parent = sum(objs)
    _effectively_at_rest(parent) &&
        return _reference_plane_orientation(reference_topology, objs, CurrentFrame())
    isobar = _sum_line_momentum(reference_topology, objs, _root_isobar_line(reference_topology))
    isobar_cm = transform_to_cmf(isobar, parent)
    return (α = azimuthal_angle(parent), cosβ = cos_theta(parent), γ = azimuthal_angle(isobar_cm))
end

function _reference_plane_orientation(
        reference_topology::DecayTopology,
        objs,
        initial_frame::AbstractInitialFrame,
    )
    angle_results = map(helicity_angle_programs(reference_topology; initial_frame)) do program
        _, result = apply_decay_instruction(program, objs)
        only(values(result))
    end
    return (α = angle_results[1].ϕ, cosβ = angle_results[1].cosθ, γ = angle_results[2].ϕ)
end

end
