# [Convention matching across amplitude frameworks](@id convention_matching)

```@meta
CurrentModule = CascadeDecays
EditURL = "../src/convention-matching.md"
```

An amplitude comparison is only meaningful after the two implementations agree on
their factorisation, particle order, helicity states, and kinematic frames. Agreement
of an unpolarised intensity at one phase-space point is not sufficient: signs may
cancel in the helicity sum, and a kinematic mismatch may accidentally look like a
constant normalisation.

This page records the conventions used by `CascadeDecays.jl` and a practical
procedure for translating models from other frameworks. It is informed by mappings
to
[`B2DxDK.jl`](https://github.com/RUB-EP1/B2DxDK.jl/blob/main/src/matching.jl),
[`Lc2ppiKSemileptonicModelLHCb.jl`](https://github.com/mmikhasenko/Lc2ppiKSemileptonicModelLHCb.jl/blob/main/src/mapping.jl),
`ThreeBodyDecays.jl`, and the models collected in
[`amplitude-serialization`](https://github.com/RUB-EP1/amplitude-serialization).
The same issues have appeared in mappings of ``\Xi_b\to pK K``,
``\Lambda_b\to pK\gamma``, and ``X\to\pi\pi\pi`` amplitudes.

The comparison target should be the **complex amplitude for every external-helicity
configuration**, evaluated at several phase-space points. A constant complex ratio
may be absorbed into a chain coupling. A ratio that depends on masses, angles, or
helicities identifies a convention mismatch and must not be hidden in a fitted
coefficient.

## The factorisation contract

For one cascade chain, `CascadeDecays.jl` evaluates the schematic amplitude

```math
\mathcal A_{\lambda_\mathrm{ext}}(x) =
\sqrt{\prod_R(2J_R+1)}
\prod_R P_R(\sigma_R;\theta_R)
\sum_{\lambda_\mathrm{int}}
\prod_v
D^{J_0*}_{\lambda_0,\lambda_1-\lambda_2}(\phi_v,\theta_v,0)
H^v_{\lambda_1\lambda_2}(x).
```

The objects in this expression have deliberately separate responsibilities:

| Object | May depend on | Must not acquire implicitly |
|---|---|---|
| Propagator ``P_R`` | Routed invariant ``\sigma_R`` and its own parameters, including fixed channel thresholds | The cascade vertex, production process, or event-dependent masses of the attached children |
| Vertex form factor ``F_v`` | The actual squared masses ``(m_0^2,m_1^2,m_2^2)`` at that vertex | The nominal mass of a propagator on line 0, 1, or 2 |
| Recoupling ``h`` | Ordered helicities, spins, and LS or helicity couplings | The particle-2 state phase applied by this package |
| Wigner matrix ``D`` | The angle measured in the local parent helicity frame | An angle reconstructed from vectors expressed in different frames |

This separation is more than software organisation. It determines which factors
are universal properties of a line and which belong to a particular occurrence of
a decay vertex.

### Propagators describe propagation

A [`Propagator`](@ref) is called as

```math
P_R(\sigma_R;\theta_R).
```

It receives only the event-dependent line invariant ``\sigma_R``. Parameters such
as pole mass, couplings, and the masses of every channel entering a self-energy or
energy-dependent width belong to ``\theta_R``. A multichannel propagator must retain
all of those thresholds even when the observed cascade contains only one of its
decay modes.

!!! danger "Do not infer a width from the attached cascade daughters"
    Replacing the propagator's channel masses by the masses of the particles at the
    current cascade vertex makes propagation depend on how the line happens to be
    used. This is a common implementation error. The two sets of masses can be
    numerically equal for a simple one-channel model, but their conceptual roles
    remain different and cease to agree for coupled channels, effective channels,
    or reuse of the same propagator in another process.

For reference, the `HadronicLineshapes.jl` Breit--Wigner used by this ecosystem has
the denominator

```math
m_R^2-\sigma-i\,m_R\Gamma(\sigma),
```

and, for one two-body channel, an energy-dependent width equivalent to

```math
\Gamma(\sigma)=\Gamma_R\,
\frac{m_R}{\sqrt{\sigma}}\,
\frac{q(\sigma)}{q_R}
\left[\frac{F_\ell(q(\sigma))}{F_\ell(q_R)}\right]^2,
\qquad
q_R=q(m_R;m_a,m_b).
```

Both the sign of the imaginary term and the precise width normalisation must be
checked when importing a lineshape. Because the ``F_\ell`` defined below already
contains the threshold power ``q^\ell``, adding a separate
``(q/q_R)^{2\ell}`` would count it twice.

### Vertex form factors describe a decay vertex

The form factor in [`Vertex`](@ref) is evaluated as

```math
F_v(m_0^2,m_1^2,m_2^2),
```

using the **actual event-dependent masses** of the ordered process
``0\to1+2``. Its breakup momentum is

```math
q(m_0;m_1,m_2)=
\frac{\sqrt{\lambda(m_0^2,m_1^2,m_2^2)}}{2m_0},
\qquad
\lambda(a,b,c)=a^2+b^2+c^2-2ab-2ac-2bc.
```

The `HadronicLineshapes.BlattWeisskopf{ℓ}(d)` functor used by this package means

```math
F_\ell(q;d)=\sqrt{\frac{z^{2\ell}}{\chi_\ell(z^2)}},
\qquad z=dq,
```

with, for example,

```math
\chi_0=1,\qquad
\chi_1=1+z^2,\qquad
\chi_2=9+3z^2+z^4,\qquad
\chi_3=225+45z^2+6z^4+z^6.
```

Thus this object includes the ``q^\ell`` threshold behaviour and tends to one at
large real ``q``. It is **not normalised at a resonance pole**: in general
``F_\ell(q_R)\ne1``. Other frameworks use the name “Blatt--Weisskopf factor” for
only the inverse-polynomial part and place ``q^\ell`` elsewhere. Compare formulas,
not class names.

A frequently used resonance-normalised convention is

```math
\widehat F_\ell(q)=\frac{F_\ell(q)}{F_\ell(q_0)},
```

where ``q_0`` is evaluated with a nominal resonance mass. This convention mixes a
pole parameter into the vertex normalisation; `CascadeDecays.jl` intentionally does
not do that. If a source chain uses ``c_\mathrm{src}\widehat F_\ell(q)`` while the
target vertex uses ``c_\mathrm{CD}F_\ell(q)``, then

```math
c_\mathrm{CD}=\frac{c_\mathrm{src}}{F_\ell(q_0)}.
```

Apply this conversion separately at every production and decay vertex. The explicit
`nominal_vertex_matching_factor` in the `B2DxDK.jl` mapping is an example of this
operation.

!!! warning "Pole normalisation may not exist"
    At a threshold, or when a nominal intermediate mass lies outside the physical
    range of its production vertex, ``q_0`` can vanish or become complex. A
    pole-normalised form factor is then singular or prescription-dependent. Do not
    silently take an absolute value, clip ``q_0``, or keep only its real part. Match
    the source framework's analytic continuation or remove the pole normalisation
    explicitly.

Keep lineshape and vertex factors distinct during this conversion. A Breit--Wigner
width may already contain a pole-normalised barrier ratio, while the observed decay
vertex contributes another amplitude-level form factor. Whether both factors are
intended is part of the model definition.

## The particle-2 state phase

At an ordered binary vertex ``0\to1+2``, `CascadeDecays.jl` distinguishes the raw
recoupling amplitude ``h`` from the helicity coupling ``H`` used in the amplitude:

```math
H_{\lambda_1\lambda_2}=
h_{\lambda_1\lambda_2}
(-1)^{j_2-\lambda_2}
F_v(m_0^2,m_1^2,m_2^2).
```

The exponent is always an integer for an allowed helicity. The phase is the
Jacob--Wick state convention for particle 2, whose momentum is opposite to particle
1. It is applied automatically by `_particle_two_phase`.

The distinction matters because parity relations and LS couplings are relations on
``h``, before this state phase is applied. In the Clebsch--Gordan convention used by
`ThreeBodyDecays.RecouplingLS`, a single ``(\ell,s)`` contribution is

```math
h_{\lambda_1\lambda_2}^{\ell s}=
\sqrt{\frac{2\ell+1}{2j_0+1}}
\langle j_1,\lambda_1;j_2,-\lambda_2
\mid s,\lambda_1-\lambda_2\rangle
\langle \ell,0;s,\lambda_1-\lambda_2
\mid j_0,\lambda_1-\lambda_2\rangle.
```

Do not impose an LS or parity relation directly on ``H`` as though the particle-2
phase were absent. Conversely, if a source framework already includes this phase in
its stored helicity couplings, do not apply it a second time.

Child order is therefore physical convention data. Swapping `(1, 2)` to `(2, 1)` is
not a harmless tree rewrite: it changes the helicity difference, the particle-2
phase, the local frame path, and generally the coupling map. Fix ordered topology
labels before translating couplings.

## Helicity frames and angle construction

Let ``R(\phi,\theta,0)=R_z(\phi)R_y(\theta)`` orient the ``z`` axis along a
particle momentum, and let ``B_z(\xi)`` be the canonical boost. The helicity-frame
map is the inverse standard transformation,

```math
L_1^{-1}(p)=[R(\phi,\theta,0)B_z(\xi)]^{-1}.
```

For particle 2 the state convention contains the extra half-turn,

```math
L_2^{-1}(p)=[R(\phi_-,\theta_-,0)R_y(\pi)B_z(\xi)]^{-1},
```

where ``(\phi_-,\theta_-)`` describe the opposite spatial momentum. To remove any
matrix-order ambiguity, the operations applied by `InstructionalDecayTrees.jl` are

| Path step | Active operations applied to every four-vector, in code order |
|---|---|
| `ToHelicityFrame` | ``R_z(-\phi)`` then ``R_y(-\theta)`` then ``B_z(-\xi)`` |
| `ToHelicityFrameParticle2` | ``R_z(-\phi_-)`` then ``R_y(-\theta_-)`` then ``R_y(-\pi)`` then ``B_z(-\xi)`` |

At every vertex, the package first follows the topology from the root into the
parent's local helicity frame. It then measures ``(\cos\theta,\phi)`` from the
momentum of ordered child 1 in that frame. A path entering child 2 uses
`ToHelicityFrameParticle2`. The vertex rotation is

```math
D^{j_0*}_{\lambda_0,\lambda_1-\lambda_2}(\phi,\theta,0).
```

Conjugation of ``D``, Euler-angle order, active versus passive rotations, azimuth
range, and the definition of child 1 must all be checked against the source.

The data flow should remain explicit:

```text
final-state four-vectors
        ↓  topology-ordered boosts and rotations
line invariants + one local angle pair per vertex
        ↓  propagators, vertex factors, and helicity contractions
complex external-helicity amplitude
```

Computing a decay-plane angle from a cross product is safe only when both vectors
are expressed in the intended common frame. Shortcuts that combine vectors from
different rest frames do not implement the helicity construction. They can agree
for two- and three-body systems because the kinematics are planar and there is only
one non-trivial decay plane. In larger cascades, successive non-collinear boosts
produce Wigner rotations and several independent planes. For four or more final
particles, stop and validate every local ``(\cos\theta,\phi)`` for every topology;
do not extrapolate a three-body cross-check.

When amplitudes with different topologies are added coherently, their external
spin states must also be transported into one reference topology. Use
[`KinematicTask`](@ref) and request the relevant `wigner_finals`; agreement of
chain-local angles alone does not validate these relative Wigner rotations.

## A reproducible matching procedure

Match from kinematics outward. At each stage, compare several generic phase-space
points, including non-zero azimuths and points away from resonance poles.

1. **Freeze topology and labels.** Record the ordered bracket tree, final-particle
   order, spin labels, external-helicity axis order, and reference topology.
2. **Compare scalar kinematics.** Check every line invariant and every vertex's
   ordered ``(m_0^2,m_1^2,m_2^2)``.
3. **Compare local frames.** Check ``\cos\theta`` and ``\phi`` one vertex at a time.
   Use four-vectors transformed by the actual frame path, not reconstructed plane
   normals.
4. **Compare each propagator.** Scan ``P_R(\sigma)`` as a complex function. Record
   denominator sign, width convention, channel masses, analytic continuation, and
   any multichannel coupling normalisation.
5. **Compare each form factor.** Write both definitions as functions of breakup
   momentum. Identify threshold powers, radius units, and every pole-normalisation
   factor.
6. **Compare one vertex.** Set propagators and all other vertices to one. Enumerate
   allowed helicities and compare ``h``, the particle-2 phase, ``H``, and the Wigner
   matrix separately.
7. **Compare one complete chain.** Test the complex amplitude for every external
   helicity. Include the package factor ``\sqrt{2J_R+1}`` for every internal line;
   some frameworks attach this factor to propagators or couplings instead.
8. **Compare topology alignment.** With more than one topology, test the external
   Wigner rotations relative to the fixed reference topology.
9. **Only then compare the coherent model.** Translate chain coefficients, choose a
   common reference coefficient and phase, and finally compare intensities.

A serialized validation block should contain enough information to reproduce this
test: a fully specified four-vector point, particle order, reference topology, and
complex amplitudes per chain and helicity. One point is useful for integrity; several
points are required to distinguish constant conversions from kinematic errors.

## Diagnosing the ratio

For a source amplitude ``A_s`` and a target amplitude ``A_t``, inspect
``r=A_t/A_s`` before changing couplings.

| Observed behaviour of ``r`` | First place to inspect |
|---|---|
| Constant for one chain and all kinematics/helicities | Coupling convention, pole form-factor normalisation, or ``\sqrt{2J+1}`` factors |
| Depends only on an invariant mass | Lineshape, running width, breakup momentum, or form factor |
| Alternating signs between particle-2 helicities | Missing or duplicated ``(-1)^{j_2-\lambda_2}`` |
| Complex phase varying with azimuth | Wigner-``D`` conjugation, Euler order, or local angle frame |
| Agreement per topology but not in the coherent sum | Reference topology or external Wigner alignment |
| Agreement only after helicity summation | State phases, helicity ordering, or accidental cancellation |
| Different constant for each resonance | Resonance-dependent pole normalisation or a model-specific rephasing |

Model-specific “magic signs” can be legitimate translations of a published
coupling convention, but they should be the **last** layer and must be documented
per chain. They are not a substitute for resolving a kinematic or helicity-dependent
mismatch.

## Related references

- M. Jacob and G. C. Wick, [*On the General Theory of Collisions for Particles
  with Spin*](https://doi.org/10.1016/0003-4916(59)90051-X), Annals of Physics 7
  (1959) 404--428.
- [Amplitude computation](@ref amplitude_computation) documents the evaluation
  pipeline and internal-line spin normalisation.
- [Routing four-vectors](@ref kinematic_tasks) documents `KinematicTask`, frame
  programs, and external Wigner alignments.
- [Topology and numbering](@ref notation) documents the ordered-child convention.
- [Cross-checking with ThreeBodyDecays](@ref cascade_vs_dpd) gives a numerical
  three-body comparison in a fixed reference topology.
