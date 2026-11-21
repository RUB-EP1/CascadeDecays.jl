# Three pions: pulling out the (12)3 orientation


This page derives the angular identity used in [the COMPASS basis
note](compass-basis.md). That note also covers reflectivity and the
executable checks. The amplitudes are `CascadeDecays` cascade
amplitudes: each decay chain is written in its own helicity frames. The
goal is to rewrite the second chain in the orientation of the first. The
rotation needed for this is the DPD parent Wigner angle, so DPD notation
from `ThreeBodyDecays` enters at exactly that step.

Label the particles $1=\pi^-$, $2=\pi^+$, $3=\pi^-$, and take the input
four-vectors in their common rest frame (the GJ frame). The physical,
ordered `CascadeDecays` chains are $\tau_{12}=((1,2),3)$ and
$\tau_{32}=((3,2),1)$. The negative pion comes first in both, as in the
COMPASS definition of the isobar decay angles. DPD admits only cyclic
orderings, so its chain for the same pair is $\tau_{23}=((2,3),1)$.
Relating $\tau_{32}$ to $\tau_{23}$ produces the sign $(-1)^S$ derived
below. All statements concern angular factors, before isobar couplings,
lineshapes, or Bose symmetrization.

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
&D^{J*}_{M\lambda}(\phi_{\rm GJ}^{32},\theta_{\rm GJ}^{32},0)
D^{S*}_{\lambda0}(\phi_{\rm TY}^{32},\theta_{\rm TY}^{32},0)\\
&\quad=\sum_{\lambda'=-J}^{J}
D^{J*}_{M\lambda'}(\Omega_{12})
d^J_{\lambda'\lambda}(\zeta)
(-1)^\lambda d^S_{\lambda0}(\theta_{\rm TY}^{32})\\
&\quad=(-1)^S\sum_{\lambda'=-J}^{J}
D^{J*}_{M\lambda'}(\Omega_{12})
d^J_{\lambda'\lambda}(\zeta)
d^S_{\lambda0}(\theta_{\rm TY}^{23}).
\end{aligned}
$$

Three details are easy to get wrong. The inner angles on the left are
those of **32**. The remaining small-$d$ has spin **$S$**. Keeping the
negative-pion-first angle $\theta^{32}$ on the right requires
**$(-1)^\lambda$**. The sum over $\lambda'$ is essential. In a full
amplitude, sum these fixed-$\lambda$ identities with the appropriate
helicity couplings.

## Why the extra phase appears

For a noncollinear three-body event,

$$
O_{12}^{-1}O_{23}=R_y(\zeta),\qquad
O_{12}^{-1}O_{32}=R_y(\zeta)R_z(\pi).
$$

The cyclic (12) and (23) frames have the same oriented plane normal.
Choosing particle 3 first in (32) reverses that normal and the
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
\theta_{\rm TY}^{23}=\pi-\theta_{\rm TY}^{32},\qquad
\phi_{\rm TY}^{23}=\phi_{\rm TY}^{32}+\pi\pmod{2\pi},
$$

and, for integer $S,\lambda$,

$$
d^S_{\lambda0}(\pi-\theta)=(-1)^{S+\lambda}d^S_{\lambda0}(\theta).
$$

Explicitly substituting $\theta_{32}=\pi-\theta_{23}$ gives

$$
\begin{aligned}
(-1)^\lambda d^S_{\lambda0}(\theta_{32})
&=(-1)^\lambda d^S_{\lambda0}(\pi-\theta_{23})\\
&=(-1)^\lambda(-1)^{S+\lambda}d^S_{\lambda0}(\theta_{23})\\
&=(-1)^{S+2\lambda}d^S_{\lambda0}(\theta_{23})\\
&=(-1)^S d^S_{\lambda0}(\theta_{23}).
\end{aligned}
$$

Here $\lambda$ is integer, so $(-1)^{2\lambda}=1$. In particular, an
odd-spin isobar changes sign when converting its ordered angular factor
from (32) to (23). This is an ordering convention; any simultaneous
conversion of couplings must use the same convention.

## The two-chain amplitude with a common orientation

Let $H^{ij}_\lambda$ denote the parent helicity coupling for an isobar
of spin $S$ and a spin-zero spectator, and $h^{ij}=H^{ij\to\pi\pi}_{00}$
the isobar decay coupling. Any scalar lineshape can multiply the
corresponding chain throughout. The couplings below retain their
original (12) and (32) ordering; the explicit $(-1)^S$ performs the
angular conversion to (23). For the two identical negative pions, add
the two chains with a plus sign:

$$
\begin{aligned}
\mathcal A_M={}&\sum_\lambda
D^{J*}_{M\lambda}(\phi_{\rm GJ}^{12},\theta_{\rm GJ}^{12},0)
H^{12}_\lambda
D^{S*}_{\lambda0}(\phi_{\rm TY}^{12},\theta_{\rm TY}^{12},0)h^{12}\\
&+\sum_\lambda
D^{J*}_{M\lambda}(\phi_{\rm GJ}^{32},\theta_{\rm GJ}^{32},0)
H^{32}_\lambda
D^{S*}_{\lambda0}(\phi_{\rm TY}^{32},\theta_{\rm TY}^{32},0)h^{32}.
\end{aligned}
$$

Absorb each inner azimuth into its parent $D$, pull the (32) orientation
into the (12) frame, and perform the explicit angle substitution above:

$$
\boxed{
\mathcal A_M=\sum_{\lambda'=-J}^{J}
D^{J*}_{M\lambda'}(\phi_{\rm GJ}^{12},\theta_{\rm GJ}^{12},\phi_{\rm TY}^{12})
\left[
H^{12}_{\lambda'}d^S_{\lambda'0}(\theta_{\rm TY}^{12})h^{12}
+(-1)^S\sum_\lambda d^J_{\lambda'\lambda}(\zeta)
H^{32}_\lambda d^S_{\lambda0}(\theta_{\rm TY}^{23})h^{32}
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
+(-1)^S d^J(\zeta)H^{32}d^S(\theta_{23})h^{32}
\right],\qquad \theta_{23}=\pi-\theta_{32}.
$$

The common orientation contains all dependence on the GJ axes. The
bracket depends on internal kinematics and couplings. For identical-pion
symmetrization the two chains use the same coupling functions, evaluated
at their respective pair invariants; an overall normalization is
conventional.

## Which zeta?

This is the **parent** rotation, `wr(1,3,0)` in `ThreeBodyDecays`, with
target spectator 1 and reference spectator 3. It is negative in this
ordering:

$$
\zeta=-\arccos\!\left(\frac{\boldsymbol p_1\cdot\boldsymbol p_3}
{|\boldsymbol p_1||\boldsymbol p_3|}\right)\in(-\pi,0).
$$

It is not the final-pion Wigner rotation `wr(1,3,1)`. The final pions
have spin zero, so their spin rotations are trivial. In invariant
notation, with $s=P^2$, $\sigma_k=(P-p_k)^2$, and the Källén function
$\Lambda$,

$$
\cos\zeta=
\frac{(s+m_1^2-\sigma_1)(s+m_3^2-\sigma_3)
-2s(\sigma_2-m_1^2-m_3^2)}
{\sqrt{\Lambda(s,m_1^2,\sigma_1)\Lambda(s,m_3^2,\sigma_3)}}.
$$

The cyclic-chain convention follows the [DPD
framework](https://arxiv.org/abs/1910.04566). Here the precise sign is
fixed independently by the Cartesian rotation identities above and
checked against the installed `ThreeBodyDecays` implementation.

## Executable cross-checks

The [COMPASS notebook](compass-basis.md) executes the checks alongside
its derivation. It generates the physical events, prints the first event
and all three sets of angles, compares the Cartesian and invariant zeta,
and tests the pair-ordering identities and reflectivity factorization.
Regenerate those results locally with:

``` sh
quarto render docs/compass-basis.qmd --to gfm
```
