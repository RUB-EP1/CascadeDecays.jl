# Isospin and kaon charge-conjugation conventions


# Goal

This note fixes a sign convention for isospin doublets and charge conjugation in amplitudes with pions and kaons, with the target application

\[
\chi_{c0},\chi_{c2} \to \pi\pi K\bar K .
\]

The initial charmonium states have

\[
I=0, \qquad C=+, \qquad J^{PC}=0^{++},2^{++} .
\]

The important point is that Clebsch--Gordan coefficients fix only the isospin recoupling signs. The charge-conjugation signs of meson pairs must be fixed from the meson currents.

Regeneration note: this tutorial is checked in as generated artifacts and is
not rendered by `docs/make.jl` in CI. After editing this notebook, run

```sh
julia --project=docs docs/render-isospin-kaon-conventions.jl
```

and commit the updated generated markdown files.

# Isospin doublets

Use the light-quark doublet

\[
q=
\begin{pmatrix}
u\\ d
\end{pmatrix},
\]

with the SU(2)-conjugate doublet convention

\[
\widehat C q
= -i\sigma_2 q^*
=
\begin{pmatrix}
-\tilde d\\
\tilde u
\end{pmatrix} .
\]

Here the tilde denotes the isospin-conjugate antiquark object. This convention is the source of the relative minus sign in antikaon doublets.

For pseudoscalar kaons,

\[
K=
\begin{pmatrix}
K^+\\ K^0
\end{pmatrix},
\qquad
K^+=u\bar s,
\qquad
K^0=d\bar s .
\]

The corresponding antikaon doublet in the same convention is

\[
\bar K=
\begin{pmatrix}
-\bar K^0\\ K^-
\end{pmatrix},
\qquad
\bar K^0=s\bar d,
\qquad
K^-=s\bar u .
\]

For vector kaons,

\[
K^*=\begin{pmatrix}
K^{*+}\\ K^{*0}
\end{pmatrix},
\qquad
\bar K^*=\begin{pmatrix}
-\bar K^{*0}\\ K^{*-}
\end{pmatrix} .
\]

These are isospin conventions. They do not yet determine the physical charge-conjugation sign of the meson state.

# Charge conjugation must be applied to currents

A kaon is not just a light-quark doublet. It is represented by a quark bilinear. Therefore the safe operation is

\[
\text{quark-field convention}
\longrightarrow
\text{meson current}
\longrightarrow
\text{meson-state phase}.
\]

For a bilinear current

\[
J_\Gamma=\bar q_1\Gamma q_2,
\]

charge conjugation gives

\[
C\,J_\Gamma\,C^{-1}
=
\eta_\Gamma\,\bar q_2\Gamma q_1,
\]

where

\[
\mathcal C^{-1}\Gamma\mathcal C
=
\eta_\Gamma\,\Gamma^T .
\]

The most useful cases are

\[
\begin{array}{c|c|c}
\text{state type} & \Gamma & \eta_\Gamma \\
\hline
0^- & i\gamma_5 & +1 \\
0^+ & 1 & +1 \\
1^- & \gamma_\rho & -1 \\
1^{++}\text{-type axial} & \gamma_\rho\gamma_5 & +1 \\
1^{+-}\text{-type axial} & \sigma_{\rho\sigma}\gamma_5\;\text{or equivalent} & -1 \\
\end{array}
\]

The sign for the axial vector current is important. It is not caused by the absence or presence of an explicit factor of \(i\). Charge conjugation is unitary, so it does not complex-conjugate \(i\). The sign follows from the Dirac matrix identity.

For the vector current,

\[
\mathcal C^{-1}\gamma_\rho\mathcal C=-\gamma_\rho^T,
\]

hence

\[
\eta_{\gamma_\rho}=-1.
\]

For the axial-vector current,

\[
\mathcal C^{-1}(\gamma_\rho\gamma_5)\mathcal C
=
(-\gamma_\rho^T)(+\gamma_5^T)
=-\gamma_\rho^T\gamma_5^T .
\]

But

\[
(\gamma_\rho\gamma_5)^T
=
\gamma_5^T\gamma_\rho^T
=-\gamma_\rho^T\gamma_5^T,
\]

because \(\{\gamma_\rho,\gamma_5\}=0\). Therefore

\[
\mathcal C^{-1}(\gamma_\rho\gamma_5)\mathcal C
=(\gamma_\rho\gamma_5)^T,
\qquad
\eta_{\gamma_\rho\gamma_5}=+1 .
\]

So, in the above isospin convention,

\[
C K = +\bar K,
\qquad
C K^* = -\bar K^* .
\]

This is the origin of the common molecular convention

\[
|K\bar K^*;C=+\rangle
=
\frac{1}{\sqrt2}
\left(
|K\bar K^*\rangle
-
|K^*\bar K\rangle
\right),
\]

and

\[
|K\bar K^*;C=-\rangle
=
\frac{1}{\sqrt2}
\left(
|K\bar K^*\rangle
+
|K^*\bar K\rangle
\right).
\]

The relative minus in the \(C=+\) molecular state is not a Clebsch sign. It is the vector-current charge-conjugation sign.

More generally, if \(X\) is an excited kaon created by \(\bar s\Gamma_X q\), define

\[
C X = \eta_X \bar X .
\]

Then a two-channel kaonic system satisfies

\[
C|K\bar X\rangle
=
\eta_X |X\bar K\rangle,
\qquad
C|X\bar K\rangle
=
\eta_X |K\bar X\rangle,
\]

up to the fixed isospin-doublet phases. The \(C\)-projected combination is therefore

\[
|K\bar X;C=\lambda\rangle
=
\frac{1}{\sqrt2}
\left(
|K\bar X\rangle
+\lambda\eta_X |X\bar K\rangle
\right),
\qquad
\lambda=\pm1 .
\]

For \(X=K^*\), \(\eta_X=-1\), hence the \(C=+\) combination has a minus sign. For an axial kaon of \(1^{++}\)-type, \(\eta_X=+1\), so the \(C=+\) combination has a plus sign. For an axial kaon of \(1^{+-}\)-type, \(\eta_X=-1\), so the \(C=+\) combination has a minus sign.

# Axial kaons and the two \(1^+\) nonets

The quark model contains two different \(1^+\) structures:

\[
{}^3P_1: \quad J^{PC}=1^{++},
\]

and

\[
{}^1P_1: \quad J^{PC}=1^{+-}.
\]

They have different charge-conjugation signs in the non-strange neutral sector,

\[
\eta_C({}^3P_1)=+1,
\qquad
\eta_C({}^1P_1)=-1.
\]

For isovector states one usually quotes \(G\)-parity,

\[
G=C\,e^{i\pi I_2}.
\]

With the usual convention, for \(I=1\),

\[
G=-C.
\]

Thus

\[
a_1(1^{++}):\quad G=-,
\qquad
b_1(1^{+-}):\quad G=+ .
\]

The strange axial kaons do not have a definite \(C\)-parity by themselves, but their signs in a \(K\bar K_1\pm K_1\bar K\) construction must be inherited from the underlying nonet/current basis. In practice:

\[
K_{1A}\sim {}^3P_1
\quad\Rightarrow\quad
\eta_{K_{1A}}=+1,
\]

whereas

\[
K_{1B}\sim {}^1P_1
\quad\Rightarrow\quad
\eta_{K_{1B}}=-1.
\]

Physical \(K_1(1270)\) and \(K_1(1400)\) are mixtures of \(K_{1A}\) and \(K_{1B}\). Therefore one should assign the charge-conjugation sign in the \((K_{1A},K_{1B})\) basis first, then rotate to the physical \(K_1\) states. Assigning a single memory-based sign directly to \(K_1(1270)\) or \(K_1(1400)\) is unsafe.

# \(K\bar K\) isospin states

Using

\[
K=
\begin{pmatrix}K^+\\K^0\end{pmatrix},
\qquad
\bar K=
\begin{pmatrix}-\bar K^0\\K^-\end{pmatrix},
\]

ordinary Clebsch--Gordan coefficients give

\[
|K\bar K;I=1,M=0\rangle
=
\frac{1}{\sqrt2}
\left(
-|K^+K^-\rangle
+|K^0\bar K^0\rangle
\right),
\]

and

\[
|K\bar K;I=0,M=0\rangle
=
-\frac{1}{\sqrt2}
\left(
|K^+K^-\rangle
+|K^0\bar K^0\rangle
\right).
\]

Thus

\[
|K^+K^-\rangle
=-\frac{1}{\sqrt2}
\left(
|I=1,M=0\rangle
+|I=0,M=0\rangle
\right),
\]

and

\[
|K^0\bar K^0\rangle
=\frac{1}{\sqrt2}
\left(
|I=1,M=0\rangle
-|I=0,M=0\rangle
\right).
\]

# Three- and four-body construction

For pions, use the isovector convention

\[
\pi^+ = |1,+1\rangle,
\qquad
\pi^0 = |1,0\rangle,
\qquad
\pi^- = |1,-1\rangle .
\]

For a three-body state such as \(\pi K\bar K\), couple first the kaon pair,

\[
\left|\left(K\bar K\right)I_{K\bar K}M_{K\bar K};\pi\right\rangle,
\]

then couple to the pion:

\[
\left|IM;I_{K\bar K}\right\rangle
=
\sum_{\mu_\pi,\mu_K,\mu_{\bar K}}
\langle \tfrac12\mu_K,\tfrac12\mu_{\bar K}|I_{K\bar K}M_{K\bar K}\rangle
\langle I_{K\bar K}M_{K\bar K},1\mu_\pi|IM\rangle
|K_{\mu_K}\bar K_{\mu_{\bar K}}\pi_{\mu_\pi}\rangle .
\]

Here \(\bar K\) means the isospin-conjugate antikaon doublet, so the physical \(K^-\) already carries the convention-dependent sign through

\[
\bar K_{-1/2}=K^- .
\]

For \(\pi\pi K\bar K\), useful couplings are

\[
(\pi\pi)_{I_{\pi\pi}}(K\bar K)_{I_{K\bar K}}
\to I=0,
\]

or sequential chains such as

\[
\pi\,(K\bar K^*)_{I,C},
\qquad
K\,(\pi\bar K)_{I},
\qquad
\bar K\,(\pi K)_{I} .
\]

The first type needs both isospin and charge-conjugation bookkeeping if the kaon system contains excited kaons.

# Constraints for \(\chi_{c0}\) and \(\chi_{c2}\)

For annihilation decays of

\[
\chi_{c0},\chi_{c2}:
\qquad
I^GJ^{PC}=0^+0^{++},\;0^+2^{++},
\]

one should enforce on the full final state

\[
I_{\rm tot}=0,
\qquad
C_{\rm tot}=+ .
\]

For topologies of the form

\[
(\pi\pi)(K\bar K),
\]

this is mostly an isospin-recouping problem. The \(K\bar K\) signs are fixed by the antikaon doublet convention above.

For topologies with excited kaons, for example

\[
K\bar K^*,
\qquad
K^*\bar K,
\qquad
K\bar K_0^*,
\qquad
K_2^*\bar K,
\]

one must additionally apply charge conjugation at the level of the corresponding meson current. The rule is

\[
C\left(\bar s\Gamma q\right)C^{-1}
=
\eta_\Gamma\,\bar q\Gamma s .
\]

The value of \(\eta_\Gamma\) depends on the spin-parity structure of the kaon resonance:

\[
\begin{array}{c|c|c}
\text{current type} & \Gamma & \eta_\Gamma \\
\hline
0^- & i\gamma_5 & +1 \\
1^- & \gamma_\rho & -1 \\
0^+ & 1 & +1 \\
1^{++}\text{-type axial} & \gamma_\rho\gamma_5 & +1 \\
1^{+-}\text{-type axial} & \sigma_{\rho\sigma}\gamma_5\;\text{or equivalent} & -1 \\
2^+ & \text{derivative tensor current} & \text{fix from chosen current}
\end{array}
\]

The last line is intentionally not compressed into a memory rule: for tensor or derivative currents, the sign should be derived from the explicit interpolating current used in the amplitude convention.

# Practical implementation recipe

1. Fix the light-quark convention:

\[
(u,d) \mapsto (-\tilde d,\tilde u).
\]

2. Define all kaon and antikaon doublets using this convention:

\[
K=(K^+,K^0)^T,
\qquad
\bar K=(-\bar K^0,K^-)^T .
\]

3. Use ordinary Clebsch--Gordan coefficients for all isospin coupling.

4. Never use Clebsches to infer the charge-conjugation sign of a \(K\bar K^*\)-type pair.

5. For each excited kaon, define the meson current \(\bar s\Gamma q\) and evaluate

\[
C(\bar s\Gamma q)C^{-1}.
\]

6. Build \(C=+\) and \(C=-\) combinations only after the current-level sign is known.

7. For \(\chi_{c0}/\chi_{c2}\to\pi\pi K\bar K\), keep only full amplitudes with

\[
I_{\rm tot}=0,
\qquad
C_{\rm tot}=+.
\]

# Minimal checks

The convention gives

\[
|K\bar K;I=0\rangle
=-\frac{1}{\sqrt2}
\left(|K^+K^-\rangle+|K^0\bar K^0\rangle\right),
\]

and

\[
|K\bar K^*;C=+\rangle
=\frac{1}{\sqrt2}
\left(|K\bar K^*\rangle-|K^*\bar K\rangle\right).
\]

These two signs have different origins:

- the first is an SU(2) doublet convention;
- the second is the vector-current sign under charge conjugation.

Keeping this separation avoids sign mistakes in \(\pi\pi K\bar K\) amplitude models.
