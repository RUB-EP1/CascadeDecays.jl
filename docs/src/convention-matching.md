# [Convention matching across amplitude frameworks](@id convention_matching)

```@meta
CurrentModule = CascadeDecays
EditURL = "../src/convention-matching.md"
```

An amplitude comparison starts by identifying how the two implementations relate
their factorisation, particle order, helicity states, and kinematic frames. Literal
agreement of all intermediate quantities is sufficient, but not necessary: two
formalisms can use different variables and still be equivalent through an explicit
mapping of amplitudes and couplings. Agreement of an unpolarised intensity at one
phase-space point is not sufficient because signs may cancel in the helicity sum and
a kinematic mismatch may accidentally look like a constant normalisation. A detailed
discussion of spin-state matching in general cascades is given by Habermann and
Mikhasenko in [*Wigner rotations for cascade
reactions*](https://inspirehep.net/literature/2827198).

This page records the conventions used by `CascadeDecays.jl` and a practical
procedure for translating models from other frameworks. It is informed by mappings
to
[`B2DxDK.jl`](https://github.com/RUB-EP1/B2DxDK.jl/blob/main/src/matching.jl),
[`Lc2ppiKSemileptonicModelLHCb.jl`](https://github.com/mmikhasenko/Lc2ppiKSemileptonicModelLHCb.jl/blob/main/src/mapping.jl),
`ThreeBodyDecays.jl`, and the models collected in
[`amplitude-serialization`](https://github.com/RUB-EP1/amplitude-serialization).
The same issues have appeared in mappings of ``\Xi_b\to pK K``,
``\Lambda_b\to pK\gamma``, and ``X\to\pi\pi\pi`` amplitudes.

The cleanest comparison target is the **complex amplitude for every
external-helicity configuration**, evaluated at several phase-space points, either
directly or after applying the proposed mapping. A constant complex ratio may be
absorbed into a chain coupling. A ratio that depends on masses, angles, or helicities
shows which part of the mapping remains incomplete.

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

The objects in this expression have deliberately separate interfaces:

| Object | Energy dependence used here | Parameters held by the object |
|---|---|---|
| Propagator ``P_R`` | Routed invariant ``\sigma_R`` through its self-energy or lineshape | Pole or bare parameters, subtraction prescription, and channel thresholds |
| Vertex form factor ``F_v`` | The actual squared masses ``(m_0^2,m_1^2,m_2^2)`` at that vertex | Its chosen functional form and scale parameters; no nominal propagator mass is supplied by the interface |
| Recoupling ``h`` | Ordered helicities and spins | LS or helicity couplings, before the package's particle-2 state phase |
| Wigner matrix ``D`` | The angle measured in the local parent helicity frame | Spin and Euler-angle convention |

This separation fixes where energy dependence enters in `CascadeDecays.jl`. Another
framework may distribute constant normalisations differently—for example between a
coupling, propagator, and pole-normalised form factor—and the matching layer should
rearrange those constants for convenience. Even energy-dependent factors can be
refactorised if the complete product is preserved. The purpose of this section is
therefore not to prescribe a unique decomposition, but to make the decomposition
used by this package explicit.

### Propagators describe propagation

A [`Propagator`](@ref) is called as

```math
P_R(\sigma_R;\theta_R).
```

It receives only the event-dependent line invariant ``\sigma_R``. Schematically, a
dressed propagator can be written as

```math
P_R(\sigma_R;\theta_R)=
\frac{1}{m_{R,\mathrm{bare}}^2-\sigma_R-\Sigma_R(\sigma_R;\theta_R)},
```

up to the chosen overall-sign and subtraction conventions. The self-energy
``\Sigma_R`` describes the virtual particle loops—or “bubbles”—that dress the state
as it propagates. Its real part shifts the mass and its imaginary part produces the
open-channel width. Pole or bare masses, channel couplings, and the masses defining
the loop thresholds belong to ``\theta_R``. A multichannel propagator can therefore
retain all of its thresholds even when the observed cascade selects only one decay
mode.

In this package, those channel masses are parameters of the propagator rather than
values routed from the attached cascade vertex. They can of course be numerically
the same in a one-channel model. A framework that constructs a lineshape together
with a particular vertex may expose the same physics through a different interface;
matching then consists of checking the resulting function of ``\sigma_R``.

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

This corresponds to retaining an absorptive self-energy
``\Sigma_R\simeq i m_R\Gamma(\sigma)`` in the denominator above. The sign of the
imaginary term and the precise width normalisation are part of the convention to be
matched. Because the ``F_\ell`` defined below already contains the threshold power
``q^\ell``, an implementation using this definition should not also introduce the
same ``(q/q_R)^{2\ell}`` dependence elsewhere.

### Vertex form factors describe a decay vertex

The form factor in [`Vertex`](@ref) is a general vertex function evaluated as

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

A form factor may represent centrifugal-barrier behaviour, a finite interaction
range, or another chosen vertex dependence. It does not have to be a
Blatt--Weisskopf function. When that choice is made, the
`HadronicLineshapes.BlattWeisskopf{L}(d)` functor used by this package means

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
    pole-normalised form factor is then singular or prescription-dependent. The
    matching should record the source prescription—analytic continuation, clipping,
    a real-valued replacement, or omission of the pole normalisation—and reproduce
    that choice explicitly.

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
1. It is applied automatically by `_particle_two_phase`. This small phase is already
present in classic presentations such as Richman's [*An Experimenter's Guide to the
Helicity Formalism*](https://inspirehep.net/literature/202987) and Martin and
Spearman's [*Elementary Particle
Theory*](https://inspirehep.net/literature/2104945), but is easy to overlook when
transcribing only the principal helicity-amplitude formulas.

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

Let ``(\phi,\theta)`` be the direction of ordered child 1 in the parent rest
frame, let ``R(\phi,\theta,0)=R_z(\phi)R_y(\theta)``, and let
``B_z(\gamma_i)`` be the canonical boost for child ``i``. The child-1 helicity-frame
map is the inverse standard transformation,

```math
L_1^{-1}=[R(\phi,\theta,0)B_z(\gamma_1)]^{-1}.
```

For child 2, the same ``R(\phi,\theta,0)`` is followed by the extra half-turn,

```math
L_2^{-1}=[R(\phi,\theta,0)R_y(\pi)B_z(\gamma_2)]^{-1}.
```

To remove any matrix-order ambiguity, the operations applied by
`InstructionalDecayTrees.jl` are

| Path step | Active operations applied to every four-vector, in code order |
|---|---|
| `ToHelicityFrame` | ``R_z(-\phi)`` then ``R_y(-\theta)`` then ``B_z(-\gamma_1)`` |
| `ToHelicityFrameParticle2` | ``R_z(-\phi)`` then ``R_y(-\theta)`` then ``R_y(-\pi)`` then ``B_z(-\gamma_2)`` |

These two standard transformations are Eq. (3) of [*Wigner rotations for cascade
reactions*](https://inspirehep.net/literature/2827198). In the notation of that
reference, the same rotation ``R(\phi,\theta,0)`` is used for both ordered daughters,
with ``R_y(\pi)`` distinguishing particle 2.

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

Cross products of momenta provide convenient definitions of decay-plane angles,
provided the frames and reference particles are stated. Those angles need not equal
the local helicity angles above for the resulting formalism to be equivalent. A
different reference particle, a reversed plane normal, or another route through
rest frames may instead lead to a mapping of the chain amplitude
``A_{\lambda_0,\lambda_1,\ldots}``, including helicity-index flips and compensating
phases in the couplings.

Appendix B of the LHCb [``\Lambda_c^+`` polarimetry
paper](https://inspirehep.net/literature/2623821) gives a concrete example: two
three-body constructions use different angle definitions, yet Eqs. (19)--(22) map
the full amplitudes through an external proton-helicity flip and explicit coupling
phases. This is the appropriate test when intermediate angles do not coincide.
For larger cascades, successive non-collinear boosts and several independent decay
planes make such equivalences less transparent. Validate every local angle or derive
the full amplitude mapping topology by topology rather than extrapolating a
three-body identity.

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
   For a direct convention match, use four-vectors transformed by the actual frame
   path. If the source uses reconstructed plane normals or other reference
   particles, test the corresponding full-amplitude mapping instead.
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
coupling convention. Document them per chain and state whether they are constant
rephasings, consequences of a helicity-index map, or part of another explicit
conversion.

## Related references

- K. Habermann and M. Mikhasenko, [*Wigner rotations for cascade
  reactions*](https://inspirehep.net/literature/2827198), Phys. Rev. D 111 (2025)
  056015.
- M. Jacob and G. C. Wick, [*On the General Theory of Collisions for Particles
  with Spin*](https://doi.org/10.1016/0003-4916(59)90051-X), Annals of Physics 7
  (1959) 404--428.
- J. D. Richman, [*An Experimenter's Guide to the Helicity
  Formalism*](https://inspirehep.net/literature/202987), CALT-68-1148 (1984).
- A. D. Martin and T. D. Spearman, [*Elementary Particle
  Theory*](https://inspirehep.net/literature/2104945), North-Holland (1970).
- LHCb collaboration, [*``\Lambda_c^+`` polarimetry using the dominant hadronic
  mode*](https://inspirehep.net/literature/2623821), JHEP 07 (2023) 228.
- [Amplitude computation](@ref amplitude_computation) documents the evaluation
  pipeline and internal-line spin normalisation.
- [Routing four-vectors](@ref kinematic_tasks) documents `KinematicTask`, frame
  programs, and external Wigner alignments.
- [Topology and numbering](@ref notation) documents the ordered-child convention.
- [Cross-checking with ThreeBodyDecays](@ref cascade_vs_dpd) gives a numerical
  three-body comparison in a fixed reference topology.
