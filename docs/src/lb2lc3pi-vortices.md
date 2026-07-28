# [Vortex manifolds in the Λb → Λc 3π model](@id lb2lc3pi_vortices)

```@meta
CurrentModule = CascadeDecays
EditURL = "../src/lb2lc3pi-vortices.md"
```

This note records the construction, numerical checks, and interactive
visualization of the zeros of the external-helicity amplitude matrix in the
toy model from [Building a full model for a decay](@ref lb2lc3pi_model). The
reproducible implementation is split between:

- [`examples/lb2lc3pi_vortices.jl`](https://github.com/RUB-EP1/CascadeDecays.jl/blob/main/examples/lb2lc3pi_vortices.jl),
  which studies the full orientation-free five-dimensional phase space;
- [`examples/lb2lc3pi_vortex_strings.jl`](https://github.com/RUB-EP1/CascadeDecays.jl/blob/main/examples/lb2lc3pi_vortex_strings.jl),
  which constructs topology-adapted coordinates and traces fixed-angle
  intersections; and
- [`examples/lb2lc3pi_vortex_lab.jl`](https://github.com/RUB-EP1/CascadeDecays.jl/blob/main/examples/lb2lc3pi_vortex_lab.jl),
  a PlutoUI and Plotly.js laboratory for changing topology and all five
  coordinates.

## What is a vortex here?

The three spin-zero pion helicities are trivial, so the model amplitude has
only the two initial and two final baryon helicities as its open indices.
Reshape it as

```math
A_{\lambda_{\Lambda_c},\lambda_{\Lambda_b}}\in\mathbb C^{2\times2}.
```

The rank-loss condition is

```math
\det A=0.
```

It is one complex equation, hence two real equations. The physical four-body
phase space has eight coordinates. Three describe an arbitrary global
orientation, and the amplitude rank is unchanged by that common rotation.
The orientation-free space therefore has dimension five. At a regular zero,

```math
\operatorname{rank}
\frac{\partial(\operatorname{Re}\det A,\operatorname{Im}\det A)}
     {\partial(u_1,\ldots,u_5)}=2,
```

so the implicit-function theorem gives a **three-dimensional zero manifold**,
not a single curve, in the full five-dimensional space.

When two relative angles are fixed, the remaining coordinates are three
invariant masses. Two real equations in that three-dimensional section
generically leave a **one-dimensional curve**. These curves are the vortex
strings drawn by the notebook. A string is therefore a fixed-angle section of
the larger three-dimensional vortex manifold.

The number of roots on one Dalitz plane is not itself a global topological
invariant. It can change when a string crosses a phase-space boundary or is
tangent to the chosen mass plane. The useful continuity tests are instead:

1. a small determinant residual;
2. rank two of the local real Jacobian normal to the zero;
3. stable continuation to neighboring slices; and
4. an integer phase winding of the determinant around a small normal loop.

The five-dimensional scan in `lb2lc3pi_vortices.jl` performs these rank,
continuation, and winding checks. The two identical negative pions exchange
some components under \(3\leftrightarrow4\).

## Topology-adapted coordinates

For any selected tree write the kinematic construction as

```text
Λb → (a,(b,c)) + d .
```

Particles \(a,b,c\) form the displayed cluster and \(d\) is the spectator.
The natural three-mass coordinates are

```math
(m_{ab},m_{ac},m_{abc}),
```

while the two relative angles are the polar angle \(\theta_a\) of particle
\(a\) in the \((abc)\) frame and the azimuth \(\phi_b\) of particle \(b\) in
the \((bc)\) frame. The notebook currently offers:

| bracket | \(a\) | \((b,c)\) | \(d\) | displayed masses | angles |
|:--|--:|:--:|--:|:--|:--|
| `((1,(3,4)),2)` | 1 | (3,4) | 2 | \(m_{13},m_{14},m_{134}\) | \(\theta_1,\phi_3\) |
| `(((1,2),3),4)` | 3 | (1,2) | 4 | \(m_{13},m_{23},m_{123}\) | \(\theta_3,\phi_1\) |
| `(((1,2),4),3)` | 4 | (1,2) | 3 | \(m_{14},m_{24},m_{124}\) | \(\theta_4,\phi_1\) |
| `((1,(2,3)),4)` | 1 | (2,3) | 4 | \(m_{12},m_{13},m_{123}\) | \(\theta_1,\phi_2\) |
| `((1,(2,4)),3)` | 1 | (2,4) | 3 | \(m_{12},m_{14},m_{124}\) | \(\theta_1,\phi_2\) |

Every topology uses the same unit-square-style internal domain. For final
masses \(m_a,m_b,m_c,m_d\),

```math
m_a+m_b+m_c < m_{abc} < m_{\Lambda_b}-m_d
```

and

```math
\begin{aligned}
m_{abc}
 &=m_a+m_b+m_c
   +u_{\rm cluster}\left(m_{\Lambda_b}-m_d-m_a-m_b-m_c\right),\\
m_{bc}
 &=m_b+m_c+t_{\rm pair}(m_{abc}-m_a-m_b-m_c),\\
\cos\chi&=2v_{\rm helicity}-1,\\
\theta_a&=\pi u_\theta,\qquad
\phi_b=\pi(2u_\phi-1).
\end{aligned}
```

Thus \(u_{\rm cluster},t_{\rm pair},v_{\rm helicity},u_\theta,u_\phi\)
all lie in \((0,1)\). At fixed \(m_{abc}\), the square
\((t_{\rm pair},v_{\rm helicity})\) maps onto its physical Dalitz domain.

The displayed masses can be calculated without reconstructing four-vectors.
In the \((bc)\) rest frame set \(s_{bc}=m_{bc}^2\) and

```math
\begin{aligned}
E_a^*&=\frac{m_{abc}^2-m_a^2-s_{bc}}{2m_{bc}},&
E_b^*&=\frac{s_{bc}+m_b^2-m_c^2}{2m_{bc}},\\
E_c^*&=\frac{s_{bc}+m_c^2-m_b^2}{2m_{bc}},&
k_a^*&=\sqrt{(E_a^*)^2-m_a^2},\\
q_b^*&=\frac{\sqrt{\lambda(s_{bc},m_b^2,m_c^2)}}{2m_{bc}}.&
\end{aligned}
```

Then

```math
\begin{aligned}
m_{ab}^2 &=m_a^2+m_b^2+
  2(E_a^*E_b^*+k_a^*q_b^*\cos\chi),\\
m_{ac}^2 &=m_a^2+m_c^2+
  2(E_a^*E_c^*-k_a^*q_b^*\cos\chi).
\end{aligned}
```

These expressions also generate the two Dalitz-border branches by taking
\(\cos\chi\to\pm1\).

## Aligned four-vector construction

The numerical event is built by nested two-body decays:

1. In the \(\Lambda_b\) rest frame, the \((abc)\) momentum is along \(+z\)
   and the spectator momentum is along \(-z\).
2. In the \((abc)\) rest frame, \(a\) lies in the \(xz\) plane at
   \(\theta_a\).
3. In the \((bc)\) rest frame, \(b\) has polar angle \(\chi\) and azimuth
   \(\phi_b\), while \(c\) is back-to-back.
4. A longitudinal boost and a \(y\) rotation place \(b,c\) in the cluster
   frame. A second longitudinal boost returns all three particles to the
   parent frame.

An exactly axial event is a coordinate pole for some helicity-frame azimuths.
Evaluating there can produce apparent determinant zeros that disappear after
an arbitrarily small common rotation. The implementation therefore applies
one fixed generic global rotation after constructing the aligned event.
Masses, relative angles, amplitude rank, and the zero set are unchanged.

For every offered topology the implementation was checked numerically to
reconstruct \(m_{ab}\), \(m_{ac}\), and \(m_{abc}\) from the four-vectors to
better than \(10^{-10}\,\mathrm{GeV}\), and to conserve the parent
four-momentum at the same scale.

## Root solving and string continuation

The code solves the scale-free function

```math
f=\frac{\det A}{\sum_{\lambda',\lambda}
  |A_{\lambda'\lambda}|^2}.
```

This removes much of the Breit–Wigner dynamic range without moving any
interior zero. On each fixed-\(m_{abc}\) plane, a damped Newton iteration
solves

```math
\operatorname{Re}f(t_{\rm pair},v_{\rm helicity})=0,\qquad
\operatorname{Im}f(t_{\rm pair},v_{\rm helicity})=0.
```

A square seed grid discovers the roots on the first plane. Roots on every
subsequent plane are used as continuation seeds, supplemented by a smaller
discovery grid. Nearest roots in the square chart are assigned the same branch
identifier only on adjacent slices; this prevents unrelated zeros from being
joined visually.

For the original `((1,(3,4)),2)` chart at
\(\theta_1=60^\circ,\phi_3=90^\circ\), a dense benchmark over
\(5.10<m_{134}<5.478\,\mathrm{GeV}\) used 150 planes and retained 518 root
intersections. The largest residual was below \(2\times10^{-9}\). A single
plane at \(m_{134}=5.35\,\mathrm{GeV}\) contains six roots, each reproducible
to approximately \(10^{-12}\) in the normalized determinant. These numbers
describe discretized intersections, not 518 separate strings.

The interactive default is deliberately lighter: seven planes in the
high-\(m_{abc}\) window. It completes in tens of seconds on a laptop and,
for the default chart and angles, returns 14 intersections with a maximum
residual of about \(1.4\times10^{-10}\). Moving an angle or changing topology
can alter the visible intersection count, so the controls are wrapped in a
confirmation button rather than recomputing at every intermediate slider
position.

## Running the Pluto laboratory

Instantiate the examples environment once:

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

Then launch the checked-in notebook:

```sh
julia --project=examples -e \
  'using Pluto; Pluto.run(notebook="examples/lb2lc3pi_vortex_lab.jl")'
```

The Plotly.js view shows several gray Dalitz borders in the three-mass volume,
red determinant zeros grouped by continuation branch, and a gold marker for
the independent five-slider phase-space point. Plotly.js is loaded in the
browser, so the plot remains rotatable and zoomable without adding a Julia
plotting backend to the notebook.

For a noninteractive dense export:

```sh
VORTEX_THETA1_DEG=60 \
VORTEX_PHI3_DEG=90 \
VORTEX_NSLICES=150 \
VORTEX_NGRID=17 \
VORTEX_MR_LO=5.10 \
VORTEX_MR_HI=5.478 \
VORTEX_JSON=/tmp/lb2lc3pi-vortices.json \
julia --project=examples examples/lb2lc3pi_vortex_strings.jl
```

The notebook has been validated both as an ordinary Julia script and by
opening and evaluating all cells through Pluto's `ServerSession`; Pluto
reported a complete reactive topology with no errored cells.
