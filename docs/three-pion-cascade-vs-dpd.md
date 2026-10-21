# Three pions: pulling out the (12)3 orientation


Label the particles $1=\pi^+$, $2=\pi^-$, $3=\pi^-$, and take the input
four-vectors in their common rest frame (the GJ frame). The ordered
chains are $\tau_{12}=((1,2),3)$ and $\tau_{13}=((1,3),2)$: the positive
pion defines the inner angles in both. DPD uses $\tau_{31}=((3,1),2)$
instead. The completed derivation, including reflectivity, is collected
in [the COMPASS basis note](compass-basis.md). All statements below
concern angular factors, before isobar couplings, lineshapes, or Bose
symmetrization.

The convention used by `PartialWaveFunctions`, and conjugated by
`CascadeDecays`, is

$$
D^{j*}_{mn}(\alpha,\beta,\gamma)
=e^{im\alpha}d^j_{mn}(\beta)e^{in\gamma},\qquad
d^j_{mn}(\beta)=\langle jm|e^{-iJ_y\beta}|jn\rangle.
$$

Write
$\Omega_{ij}=(\phi_{\rm GJ}^{ij},\theta_{\rm GJ}^{ij},\phi_{\rm TY}^{ij})$
and
$O_{ij}=R_z(\phi_{\rm GJ}^{ij})R_y(\theta_{\rm GJ}^{ij})R_z(\phi_{\rm TY}^{ij})$.
The $z$ axis points along the pair momentum, and the first particle of
the ordered pair has positive $x$ in the resulting decay-plane frame.

## Exact identity

For fixed $\lambda$ (with $|\lambda|\leq\min(J,S)$),

$$
\begin{aligned}
&D^{J*}_{M\lambda}(\phi_{\rm GJ}^{13},\theta_{\rm GJ}^{13},0)
D^{S*}_{\lambda0}(\phi_{\rm TY}^{13},\theta_{\rm TY}^{13},0)\\
&\quad=\sum_{\lambda'=-J}^{J}
D^{J*}_{M\lambda'}(\Omega_{12})
d^J_{\lambda'\lambda}(\zeta)
(-1)^\lambda d^S_{\lambda0}(\theta_{\rm TY}^{13})\\
&\quad=(-1)^S\sum_{\lambda'=-J}^{J}
D^{J*}_{M\lambda'}(\Omega_{12})
d^J_{\lambda'\lambda}(\zeta)
d^S_{\lambda0}(\theta_{\rm TY}^{31}).
\end{aligned}
$$

There are three corrections to the proposed expression: the inner angles
on the left must be those of **13**, the remaining small-$d$ has spin
**$S$**, and using the positive-pion-first angle on the right requires
**$(-1)^\lambda$**. The sum over $\lambda'$ is essential. In a full
amplitude, sum these fixed-$\lambda$ identities with the appropriate
helicity couplings.

## Why the extra phase appears

For a noncollinear three-body event,

$$
O_{12}^{-1}O_{31}=R_y(\zeta),\qquad
O_{12}^{-1}O_{13}=R_y(\zeta)R_z(\pi).
$$

The cyclic (12) and (31) frames have the same oriented plane normal.
Choosing particle 1 first in (13) reverses that normal and the
transverse $x$ axis. The last $R_z(\pi)$ contributes
$e^{i\lambda\pi}=(-1)^\lambda$. Meanwhile the inner azimuth is absorbed
into the last Euler angle of the parent matrix:

$$
D^{J*}_{M\lambda}(\phi,\theta,0)
D^{S*}_{\lambda0}(\varphi,\vartheta,0)
=D^{J*}_{M\lambda}(\phi,\theta,\varphi)d^S_{\lambda0}(\vartheta).
$$

Reversing the pair gives

$$
\theta_{\rm TY}^{31}=\pi-\theta_{\rm TY}^{13},\qquad
\phi_{\rm TY}^{31}=\phi_{\rm TY}^{13}+\pi\pmod{2\pi},
$$

and, for integer $S,\lambda$,

$$
d^S_{\lambda0}(\pi-\theta)=(-1)^{S+\lambda}d^S_{\lambda0}(\theta).
$$

Explicitly substituting $\theta_{13}=\pi-\theta_{31}$ gives

$$
\begin{aligned}
(-1)^\lambda d^S_{\lambda0}(\theta_{13})
&=(-1)^\lambda d^S_{\lambda0}(\pi-\theta_{31})\\
&=(-1)^\lambda(-1)^{S+\lambda}d^S_{\lambda0}(\theta_{31})\\
&=(-1)^{S+2\lambda}d^S_{\lambda0}(\theta_{31})\\
&=(-1)^S d^S_{\lambda0}(\theta_{31}).
\end{aligned}
$$

Here $\lambda$ is integer, so $(-1)^{2\lambda}=1$. In particular, an
odd-spin isobar changes sign when converting its ordered angular factor
from (13) to (31). This is an ordering convention; any simultaneous
conversion of couplings must use the same convention.

## The two-chain amplitude with a common orientation

Let $H^{ij}_\lambda$ denote the parent helicity coupling for an isobar
of spin $S$ and a spin-zero spectator, and $h^{ij}=H^{ij\to\pi\pi}_{00}$
the isobar decay coupling. Any scalar lineshape can multiply the
corresponding chain throughout. The couplings below retain their
original (12) and (13) ordering; the explicit $(-1)^S$ performs the
angular conversion to (31). For the two identical negative pions, add
the two chains with a plus sign:

$$
\begin{aligned}
\mathcal A_M={}&\sum_\lambda
D^{J*}_{M\lambda}(\phi_{\rm GJ}^{12},\theta_{\rm GJ}^{12},0)
H^{12}_\lambda
D^{S*}_{\lambda0}(\phi_{\rm TY}^{12},\theta_{\rm TY}^{12},0)h^{12}\\
&+\sum_\lambda
D^{J*}_{M\lambda}(\phi_{\rm GJ}^{13},\theta_{\rm GJ}^{13},0)
H^{13}_\lambda
D^{S*}_{\lambda0}(\phi_{\rm TY}^{13},\theta_{\rm TY}^{13},0)h^{13}.
\end{aligned}
$$

Absorb each inner azimuth into its parent $D$, pull the (13) orientation
into the (12) frame, and perform the explicit angle substitution above:

$$
\boxed{
\mathcal A_M=\sum_{\lambda'=-J}^{J}
D^{J*}_{M\lambda'}(\phi_{\rm GJ}^{12},\theta_{\rm GJ}^{12},\phi_{\rm TY}^{12})
\left[
H^{12}_{\lambda'}d^S_{\lambda'0}(\theta_{\rm TY}^{12})h^{12}
+(-1)^S\sum_\lambda d^J_{\lambda'\lambda}(\zeta)
H^{13}_\lambda d^S_{\lambda0}(\theta_{\rm TY}^{31})h^{13}
\right].
}
$$

The internal sum runs over $|\lambda|\leq\min(J,S)$; the first term is
zero for $|\lambda'|>S$. With helicity contractions understood, this is

$$
\mathcal A=
D^{J*}(\Omega_{12})
\left[
H^{12}d^S(\theta_{12})h^{12}
+(-1)^S d^J(\zeta)H^{13}d^S(\theta_{31})h^{13}
\right],\qquad \theta_{31}=\pi-\theta_{13}.
$$

The common orientation contains all dependence on the GJ axes. The
bracket depends on internal kinematics and couplings. For identical-pion
symmetrization the two chains use the same coupling functions, evaluated
at their respective pair invariants; an overall normalization is
conventional.

## Which zeta?

This is the **parent** rotation, `wr(2,3,0)` in `ThreeBodyDecays`, with
target spectator 2 and reference spectator 3. It is positive in this
ordering:

$$
\zeta=\arccos\!\left(\frac{\boldsymbol p_2\cdot\boldsymbol p_3}
{|\boldsymbol p_2||\boldsymbol p_3|}\right)\in(0,\pi).
$$

It is not the final-pion Wigner rotation `wr(2,3,1)`. The final pions
have spin zero, so their spin rotations are trivial. In invariant
notation, with $s=P^2$, $\sigma_k=(P-p_k)^2$, and the Källén function
$\Lambda$,

$$
\cos\zeta=
\frac{(s+m_2^2-\sigma_2)(s+m_3^2-\sigma_3)
-2s(\sigma_1-m_2^2-m_3^2)}
{\sqrt{\Lambda(s,m_2^2,\sigma_2)\Lambda(s,m_3^2,\sigma_3)}}.
$$

The cyclic-chain convention follows the [DPD
framework](https://arxiv.org/abs/1910.04566). Here the precise sign is
fixed independently by the Cartesian rotation identities above and
checked against the installed `ThreeBodyDecays` implementation.

## Reproducible numerical validation

The script generates 100 random physical events at $\sqrt{s}=1.4$ GeV,
$m_\pi=0.13957039$ GeV, with seed 20261021. It draws two spatial
momenta, sets the third to minus their sum, and rescales to the chosen
total energy. This is a kinematic stress test, not a uniform phase-space
sample.

Angles come from `CascadeDecays.helicity_angle_programs` with
`CurrentFrame()`. They are checked independently using Cartesian axes
and a Lorentz boost. The test compares DPD’s invariant zeta with the
signed Cartesian rotation, then checks both exact identities for every
$J=0,\ldots,4$, $S=0,\ldots,3$, every $M$, and every allowed $\lambda$.
Collinear endpoints are excluded because their decay-plane azimuths are
undefined.

Run from the repository root:

``` sh
julia --project=. docs/validate-three-pion-dpd.jl
```

``` julia
include(joinpath(@__DIR__, "validate-three-pion-dpd.jl"))
```

    Test Summary:                    | Pass  Total  Time
    Three-pion angular factorization |  511    511  1.2s
    Seed = 20261021; events = 100; J = 0:4; S = 0:3
    First event four-vectors (E, px, py, pz), GeV:
    p1 = ( 0.511082386798,  0.113961653220,  0.455775068512, -0.144938403634)
    p2 = ( 0.374615321371, -0.126905897469,  0.005677076243,  0.323603784332)
    p3 = ( 0.514302291831,  0.012944244249, -0.461452144755, -0.178665380698)
    chain     phi_GJ         theta_GJ       phi_TY         theta_TY (radians)
    (12)3   1.598840081439  1.201522051165 -0.401452404728  1.268141169636
    (13)2  -0.044704729637  2.767519388212  1.644623660411  1.580472501202
    (31)2  -0.044704729637  2.767519388212 -1.496968993179  1.561120152388
    zeta = 1.939865882161 rad = 111.146127869224 deg
    max angles          7.105427e-15
    max missing_phase   1.365251e+00
    max mixed_TY        1.973362e+00
    max pull13          3.740991e-14
    max pull31          3.740991e-14
    max reflectivity    2.140804e-14
    max reflectivity_M0 9.711508e-15
    max reversal        9.756085e-15
    max rotation13      3.249403e-14
    max rotation31      3.339097e-14
    max zeta            2.309264e-14

    (first_event = (ps = ([0.11396165322047952, 0.45577506851214106, -0.14493840363433785, 0.5110823867976806], [-0.12690589746907208, 0.00567707624299945, 0.323603784332115, 0.3746153213710266], [0.012944244248592555, -0.4614521447551405, -0.17866538069777715, 0.5143022918312926]), a = ((phi = 1.598840081439314, theta = 1.2015220511651594), (phi = -0.4014524047283381, theta = 1.2681411696361962)), b = ((phi = -0.044704729637404114, theta = 2.7675193882118228), (phi = 1.6446236604112485, theta = 1.5804725012020748)), c = ((phi = -0.044704729637404114, theta = 2.7675193882118228), (phi = -1.4969689931785448, theta = 1.5611201523877176)), zeta = 1.9398658821605894, q = [-0.3607477965838428 2.7755575615628914e-16 0.9326634051252909; 5.204170427930421e-18 1.0000000000000002 -4.440892098500626e-16; -0.9326634051252909 -1.6653345369377348e-16 -0.3607477965838428]), errors = Dict(:reflectivity => 2.1408039610327508e-14, :pull31 => 3.74099060171673e-14, :zeta => 2.3092638912203256e-14, :reflectivity_M0 => 9.711508340058e-15, :mixed_TY => 1.97336153911434, :reversal => 9.756084828893563e-15, :rotation13 => 3.249402964114936e-14, :rotation31 => 3.3390972979997585e-14, :pull13 => 3.74099060171673e-14, :missing_phase => 1.3652508509331431…))

The validated run, including the reflectivity extension, passed all 511
assertions. Maximum absolute residuals were $3.75\times10^{-14}$ for
either pulled-out identity and $2.31\times10^{-14}$ rad for zeta.
Omitting the phase gives a maximum residual of 1.365; mixing chain-12 TY
angles into chain 13 gives 1.973.

The first event is $(E,p_x,p_y,p_z)$ in GeV:

| Particle |            $E$ |           $p_x$ |           $p_y$ |           $p_z$ |
|----------|---------------:|----------------:|----------------:|----------------:|
| 1        | 0.511082386798 |  0.113961653220 |  0.455775068512 | -0.144938403634 |
| 2        | 0.374615321371 | -0.126905897469 |  0.005677076243 |  0.323603784332 |
| 3        | 0.514302291831 |  0.012944244249 | -0.461452144755 | -0.178665380698 |

| Chain | $\phi_{\rm GJ}$ | $\theta_{\rm GJ}$ | $\phi_{\rm TY}$ | $\theta_{\rm TY}$ |
|-------|----------------:|------------------:|----------------:|------------------:|
| (12)3 |  1.598840081439 |    1.201522051165 | -0.401452404728 |    1.268141169636 |
| (13)2 | -0.044704729637 |    2.767519388212 |  1.644623660411 |    1.580472501202 |
| (31)2 | -0.044704729637 |    2.767519388212 | -1.496968993179 |    1.561120152388 |

All angles are in radians, and $\zeta=1.939865882161$ rad
$=111.146127869224^\circ$.
