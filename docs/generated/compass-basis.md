# COMPASS: three-pion meson amplitudes


This comparison with `ThreeBodyDecays.jl` specializes the [general
relations and signs](cascade-vs-dpd.md) to a meson decaying into three
spin-zero pions. It derives the common orientation, the cyclic-pair
phase, and the COMPASS reflectivity basis. The [detailed angular
cross-check](three-pion-cascade-vs-dpd.md) supplies an explicit
four-vector configuration and the angles of both topologies.

For $1=\pi^+$, $2=\pi^-$, and $3=\pi^-$, the two-chain amplitude in the
COMPASS reflectivity basis has the factorized form

``` math
\boxed{
{}^\epsilon\mathcal A_M
=\sum_{\nu=-J}^{J}{}^\epsilon D^{J*}_{M\nu}
(\phi_{\rm GJ}^{12},\theta_{\rm GJ}^{12},\phi_{\rm TY}^{12})
\mathcal O^\nu.
}
```

The reflectivity transformation acts on the parent spin projection $M$.
It leaves the internal amplitude $\mathcal O^\nu$ unchanged. The first
and third Euler angles are distinct: they are the GJ azimuth and TY
azimuth, respectively. This note fixes the ordering and phase
conventions needed to obtain that result from cascade amplitudes.

## Angles, spins, and ordering

The input momenta are in their common rest frame, with the GJ axes fixed
by production. The pair momentum defines the outer helicity axis. The
first particle of the ordered pair defines the inner TY angles in the
pair rest frame. Thus $(12)3$ and $(13)2$ both use the positive pion
first; DPD’s cyclic second chain is $(31)2$.

We use the convention implemented by `PartialWaveFunctions`:

``` math
D^{j*}_{mn}(\alpha,\beta,\gamma)
=e^{im\alpha}d^j_{mn}(\beta)e^{in\gamma},\qquad
d^j_{mn}(\beta)=\langle jm|e^{-iJ_y\beta}|jn\rangle.
```

Here $J$ is the parent spin, $S$ is the isobar spin, and the pions have
spin zero. Write $\theta_{ij}=\theta_{\rm TY}^{ij}$ and
$\Omega_{ij}=(\phi_{\rm GJ}^{ij},\theta_{\rm GJ}^{ij},\phi_{\rm TY}^{ij})$.
Let $H^{ij}_\lambda$ be the parent helicity coupling and
$h^{ij}=H^{ij\to\pi\pi}_{00}$ the isobar decay coupling. A scalar
lineshape may multiply each chain throughout the following equations.

## Pulling out the common (12) orientation

The Bose-symmetrized angular amplitude starts with a plus between the
positive-pion-first chains:

``` math
\begin{aligned}
\mathcal A_M={}&\sum_\lambda
D^{J*}_{M\lambda}(\phi_{\rm GJ}^{12},\theta_{\rm GJ}^{12},0)
H^{12}_\lambda
D^{S*}_{\lambda0}(\phi_{\rm TY}^{12},\theta_{12},0)h^{12}\\
&+\sum_\lambda
D^{J*}_{M\lambda}(\phi_{\rm GJ}^{13},\theta_{\rm GJ}^{13},0)
H^{13}_\lambda
D^{S*}_{\lambda0}(\phi_{\rm TY}^{13},\theta_{13},0)h^{13}.
\end{aligned}
```

An overall Bose normalization is conventional. Physical identical-pion
chains use the same coupling functions at the respective pair
invariants. Absorb the inner azimuth into the parent matrix using

``` math
D^{J*}_{M\lambda}(\phi,\theta,0)
D^{S*}_{\lambda0}(\varphi,\vartheta,0)
=D^{J*}_{M\lambda}(\phi,\theta,\varphi)d^S_{\lambda0}(\vartheta).
```

For
$O_{ij}=R_z(\phi_{\rm GJ}^{ij})R_y(\theta_{\rm GJ}^{ij})R_z(\phi_{\rm TY}^{ij})$,
the Cartesian rotations obey

``` math
O_{12}^{-1}O_{31}=R_y(\zeta),\qquad
O_{12}^{-1}O_{13}=R_y(\zeta)R_z(\pi).
```

Consequently,

``` math
D^{J*}_{M\lambda}(\Omega_{13})
=\sum_\nu D^{J*}_{M\nu}(\Omega_{12})
d^J_{\nu\lambda}(\zeta)(-1)^\lambda.
```

The angle is the parent rotation `ThreeBodyDecays.wr(2,3,0)`, with
reference spectator 3 and target spectator 2:

``` math
\zeta=\arccos\frac{\boldsymbol p_2\cdot\boldsymbol p_3}
{|\boldsymbol p_2||\boldsymbol p_3|}\in(0,\pi).
```

The ordered pair swap gives $\theta_{31}=\pi-\theta_{13}$ and
$\phi_{\rm TY}^{31}=\phi_{\rm TY}^{13}+\pi\pmod{2\pi}$. The small-$d$
substitution, explicitly, is

``` math
\begin{aligned}
(-1)^\lambda d^S_{\lambda0}(\theta_{13})
&=(-1)^\lambda d^S_{\lambda0}(\pi-\theta_{31})\\
&=(-1)^\lambda(-1)^{S+\lambda}d^S_{\lambda0}(\theta_{31})\\
&=(-1)^S d^S_{\lambda0}(\theta_{31}),
\end{aligned}
```

since $\lambda$ is integer. Therefore

``` math
\mathcal A_M=\sum_\nu D^{J*}_{M\nu}(\Omega_{12})\mathcal O^\nu,
```

``` math
\boxed{
\mathcal O^\nu=
H^{12}_\nu d^S_{\nu0}(\theta_{12})h^{12}
+(-1)^S\sum_\lambda d^J_{\nu\lambda}(\zeta)
H^{13}_\lambda d^S_{\lambda0}(\theta_{31})h^{13}.
}
```

The sum over $\nu$ runs from $-J$ to $J$; the first term vanishes for
$|\nu|>S$. The internal sum has $|\lambda|\leq\min(J,S)$. The couplings
retain their original (13) ordering: the explicit $(-1)^S$ already
accounts for the angular conversion to (31). The bracket depends on
internal kinematics and couplings, not the orientation of the GJ axes.

## Reflectivity convention

For parent parity $P=\pm1$, use the COMPASS convention in [Mikhasenko,
*Three-pion dynamics at COMPASS*, Eqs.
(3.10)–(3.11)](https://wwwcompass.cern.ch/compass/publications/theses/2019_phd_mikhasenko.pdf):

``` math
{}^\epsilon D^J_{M\nu}(\Omega)
=N_M\left[D^J_{M\nu}(\Omega)
-\epsilon P(-1)^{J-M}D^J_{-M,\nu}(\Omega)\right],
\qquad \epsilon=\pm1,
```

``` math
N_M=\begin{cases}
1/\sqrt2,&M>0,\\
1/2,&M=0,\\
0,&M<0.
\end{cases}
```

The coefficients are real, so conjugation gives the same combination of
$D^{J*}$. In explicit Euler angles,

``` math
{}^\epsilon D^{J*}_{M\nu}(\phi,\theta,\psi)
=N_M e^{i\nu\psi}
\left[e^{iM\phi}d^J_{M\nu}(\theta)
-\epsilon P(-1)^{J-M}e^{-iM\phi}d^J_{-M,\nu}(\theta)\right].
```

Thus the transformation combines $M$ and $-M$, including their small-$d$
rows, while retaining the same $\nu$ and third-angle phase. It is a
row-basis transformation, not just replacement of the first-angle
exponential by a sine or cosine for general $\nu$.

At $M=0$ the definition reduces to

``` math
{}^\epsilon D^{J*}_{0\nu}
=\frac{1-\epsilon P(-1)^J}{2}D^{J*}_{0\nu}.
```

Only $\epsilon=-P(-1)^J$ survives in this sign convention; the other
combination vanishes. For example, $J^P=1^+$ permits $M^\epsilon=0^+$.

## Reflectivity acts outside the internal amplitude

Apply the same linear combination to the complete amplitude:

``` math
\begin{aligned}
{}^\epsilon\mathcal A_M
&=N_M\left[\mathcal A_M-\epsilon P(-1)^{J-M}\mathcal A_{-M}\right]\\
&=\sum_\nu N_M\left[D^{J*}_{M\nu}(\Omega_{12})
-\epsilon P(-1)^{J-M}D^{J*}_{-M,\nu}(\Omega_{12})\right]\mathcal O^\nu\\
&=\boxed{\sum_\nu{}^\epsilon D^{J*}_{M\nu}(\Omega_{12})\mathcal O^\nu}.
\end{aligned}
```

In matrix notation, if $U_\epsilon$ denotes the row transformation,

``` math
U_\epsilon[D^{J*}(\Omega_{12})\mathcal O]
=[U_\epsilon D^{J*}(\Omega_{12})]\mathcal O.
```

There is no additional reflectivity factor in $d^J(\zeta)$ or
$\mathcal O$. The pair-ordering factor $(-1)^S$ remains exactly where it
was. This argument applies to each fixed $J,P$ wave; any physical parity
constraints on the helicity couplings are separate from the linear
change of basis.

## Integrating out the orientation

Hold the internal invariants fixed and use the normalized Haar measure

``` math
d\mu(\Omega_{12})=
\frac{d\Omega_{\rm GJ}^{12}\,d\phi_{\rm TY}^{12}}{8\pi^2}
=\frac{d\phi_{\rm GJ}^{12}\,d(\cos\theta_{\rm GJ}^{12})\,
d\phi_{\rm TY}^{12}}{8\pi^2}.
```

Both azimuths cover a full $2\pi$ interval and $\cos\theta$ runs from
$-1$ to $1$. This integral removes the orientation only; the internal TY
polar angles and $\zeta$ remain in $\mathcal O$ as functions of the
Dalitz invariants. There is no acceptance weighting.

Wigner orthogonality gives

``` math
\int d\mu(\Omega)\,D^{J*}_{m\lambda}(\Omega)
D^{J'}_{m'\lambda'}(\Omega)
=\frac{\delta_{JJ'}\delta_{mm'}\delta_{\lambda\lambda'}}{2J+1}.
```

At a fixed parent parity $P$, the allowed reflectivity rows form an
orthonormal change of basis. For $M>0$, their row norm is
$N_M^2(1+1)=1$. For $M=0$, the allowed sector has coefficient
$N_0[1-\epsilon P(-1)^J]=(1/2)\times2=1$; the forbidden sector is the
zero vector and is excluded from the basis. Consequently,

``` math
\int d\mu(\Omega)\,
{}^\epsilon D^{J*}_{M\lambda}(\Omega)
{}^{\epsilon'}D^{J'}_{M'\lambda'}(\Omega)
=\frac{\delta_{JJ'}\delta_{MM'}\delta_{\epsilon\epsilon'}
\delta_{\lambda\lambda'}}{2J+1}.
```

Let $a,b$ label two internal waves or models, with

``` math
\Psi^a_{JM\epsilon}
=\sum_\lambda{}^\epsilon D^{J*}_{M\lambda}(\Omega_{12})
\mathcal O_{J,a}^{\lambda}.
```

Their orientation-integrated overlap is therefore

``` math
\begin{aligned}
\int d\mu\,\Psi^a_{JM\epsilon}
\bigl(\Psi^b_{J'M'\epsilon'}\bigr)^*
&=\sum_{\lambda,\lambda'}\mathcal O_{J,a}^{\lambda}
\bigl(\mathcal O_{J',b}^{\lambda'}\bigr)^*
\int d\mu\,{}^\epsilon D^{J*}_{M\lambda}
{}^{\epsilon'}D^{J'}_{M'\lambda'}\\
&=\boxed{\frac{\delta_{JJ'}\delta_{MM'}\delta_{\epsilon\epsilon'}}{2J+1}
\sum_{\lambda=-J}^{J}\mathcal O_{J,a}^{\lambda}
\bigl(\mathcal O_{J,b}^{\lambda}\bigr)^*}.
\end{aligned}
```

For equal reflectivity, simply drop $\delta_{\epsilon\epsilon'}$. For a
single wave this becomes

``` math
\int d\mu\,|\Psi_{JM\epsilon}|^2
=\frac{1}{2J+1}\sum_{\lambda=-J}^{J}|\mathcal O_J^\lambda|^2.
```

Thus the reflectivity basis **preserves** the usual Wigner-$D$
normalization. If the wave definition includes $\sqrt{2J+1}$, namely
$\widehat\Psi=\sqrt{2J+1}\,\Psi$, its overlap instead equals the same
Kronecker deltas times $\sum_\lambda\mathcal O_a^\lambda
\mathcal O_b^{\lambda *}$, without the denominator. Any such factor
already included in the couplings must be counted only once.

These expressions assume the two waves use the same parent parity. If
their parities are $P$ and $P'$, the row-overlap factor is
$\delta_{\epsilon P,\epsilon'P'}$ instead of
$\delta_{\epsilon\epsilon'}$, for allowed rows. An orientation integral
alone does not imply an additional $\delta_{PP'}$.

## Numerical validation

Run from the repository root:

``` sh
julia --project=. docs/validate-three-pion-dpd.jl
```

The validation uses 100 random physical rest-frame events with
$m_{3\pi}=1.4$ GeV, $m_\pi=0.13957039$ GeV, and seed 20261021. It checks
the measured angles against an independent Cartesian boost, the parent
zeta against DPD invariants, both pair-ordering identities, and the full
two-chain reflectivity factorization with complex couplings. The last
check uses arbitrary couplings to test the algebra independently of
dynamical parity constraints.

For $J=0,\ldots,4$, $S=0,\ldots,3$, all allowed projections, $P=\pm1$,
and $\epsilon=\pm1$, all 511 assertions pass:

| Check                                     | Maximum absolute residual |
|:------------------------------------------|--------------------------:|
| Parent zeta, Cartesian vs invariant (rad) |      $2.31\times10^{-14}$ |
| Ordered-chain angular factorization       |      $3.75\times10^{-14}$ |
| Full amplitude in the reflectivity basis  |      $2.15\times10^{-14}$ |
| Surviving and vanishing $M=0$ sectors     |      $9.72\times10^{-15}$ |

The first event has $\zeta=1.939865882161$ rad. Omitting the ordering
phase gives residuals as large as 1.365. Noncollinear events are used,
since decay-plane azimuths are undefined at collinear endpoints.

The following output is generated when this note is rendered locally:

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
