```@meta
CurrentModule = CascadeDecays
EditURL = "../four-pion-model.qmd"
```

# [Four-pion model-building catalogue](@id four_pion_model)


This note builds a readable catalogue for

``` math
\chi_{cJ} \to \pi^+_1 \pi^-_2 \pi^+_3 \pi^-_4,
\qquad J = 0,2.
```

The first half is purely kinematical: it identifies which topology
copies are needed before any resonance model is attached. The second
half adds dynamics: toy Breit-Wigner lineshapes, LS-coupling expansion,
a `KinematicTask`, RAMBO phase-space events, and small timing checks.

The two parent hypotheses are

``` math
\chi_{c0}(0^{++}), \qquad \chi_{c2}(2^{++}).
```

The topology catalogue is common to both. The LS lists differ because
the production vertex knows the parent spin and parity.

Regeneration note: this tutorial is intentionally checked in as
generated artifacts and is not rendered by `docs/make.jl` in CI. After
editing this notebook, run

``` sh
julia --project=docs docs/render-generated-docs.jl
```

and commit the updated `docs/generated/four-pion-model.md` file.

## Setup

| parent | jp  | mass   |
|--------|-----|--------|
| chi_c0 | 0+  | 3.4147 |
| chi_c2 | 2+  | 3.5562 |

| name | jp  | mass  | width |
|------|-----|-------|-------|
| f0   | 0+  | 0.98  | 0.07  |
| rho  | 1-  | 0.775 | 0.149 |
| f2   | 2+  | 1.275 | 0.185 |

| name | jp  | mass | width |
|------|-----|------|-------|
| a1   | 1+  | 1.23 | 0.42  |
| a2   | 2+  | 1.32 | 0.11  |
| pi0  | 0-  | 1.3  | 0.25  |
| pi1  | 1-  | 1.6  | 0.3   |
| pi2  | 2-  | 1.67 | 0.26  |

## Kinematic Topology Catalogue

The cascade class starts from one positive seed,

``` math
(((1,2),3),4)
\equiv ((\pi^+_1\pi^-_2)\pi^+_3)\pi^-_4 .
```

Same-charge permutations generate the positive copies. Charge
conjugation gives one negative seed, and applying the same
charge-preserving permutations to that seed generates the negative
copies. The order of children is preserved in this class, so the
conjugate of `(1,2)` is `(2,1)`.

| stage | R3_charge | R3_pions | R2_pair | bachelor | spectator | bracket |
|----|----|----|----|----|----|----|
| positive: identity | \+ | pi1+ pi2- pi3+ | pi1+ pi2- | pi3+ | pi4- | `(((1, 2), 3), 4)` |
| positive: swap pi+ | \+ | pi3+ pi2- pi1+ | pi3+ pi2- | pi1+ | pi4- | `(((3, 2), 1), 4)` |
| positive: swap pi- | \+ | pi1+ pi4- pi3+ | pi1+ pi4- | pi3+ | pi2- | `(((1, 4), 3), 2)` |
| positive: swap pi+ and pi- | \+ | pi3+ pi4- pi1+ | pi3+ pi4- | pi1+ | pi2- | `(((3, 4), 1), 2)` |
| conjugate: identity | \- | pi2- pi1+ pi4- | pi2- pi1+ | pi4- | pi3+ | `(((2, 1), 4), 3)` |
| conjugate: swap pi+ | \- | pi2- pi3+ pi4- | pi2- pi3+ | pi4- | pi1+ | `(((2, 3), 4), 1)` |
| conjugate: swap pi- | \- | pi4- pi1+ pi2- | pi4- pi1+ | pi2- | pi3+ | `(((4, 1), 2), 3)` |
| conjugate: swap pi+ and pi- | \- | pi4- pi3+ pi2- | pi4- pi3+ | pi2- | pi1+ | `(((4, 3), 2), 1)` |

The second topology class has two two-body isobars,

``` math
\chi_{cJ} \to R_2(\pi^+\pi^-)\,R_2(\pi^+\pi^-).
```

Here reversed child order inside both two-body pairs is a bookkeeping
duplicate. This is the useful cross-check:

``` math
((1,2),(3,4)) \sim ((2,1),(4,3)),
\qquad
A_{(12)(34)} + A_{(21)(43)} = 2 A_{(12)(34)} .
```

| R2a_pair  | R2b_pair  | bracket            |
|-----------|-----------|--------------------|
| pi1+ pi2- | pi3+ pi4- | `((1, 2), (3, 4))` |
| pi1+ pi4- | pi3+ pi2- | `((1, 4), (3, 2))` |

| stage                  | raw                | canonical          |
|------------------------|--------------------|--------------------|
| Bose: identity         | `((1, 2), (3, 4))` | `((1, 2), (3, 4))` |
| Bose: swap pi+         | `((3, 2), (1, 4))` | `((1, 4), (3, 2))` |
| Bose: swap pi-         | `((1, 4), (3, 2))` | `((1, 4), (3, 2))` |
| Bose: swap pi+ and pi- | `((3, 4), (1, 2))` | `((1, 2), (3, 4))` |
| conjugate seed         | `((2, 1), (4, 3))` | `((1, 2), (3, 4))` |

## LS Catalogue

The following tables attach resonance quantum numbers and ask
`CascadeDecays.jl` for the allowed LS couplings. The construction is
still kinematic plus spin-parity bookkeeping; no lineshape parameters
enter yet.

| parent | r3 | r2 | parent -\> R3 pi | R3 -\> R2 pi | R2 -\> pi pi | n_ls | explicit_chains |
|----|----|----|----|----|----|----|----|
| chi_c0 | a1 | f0 | (L=1, S=1) | (L=1, S=0) | (L=0, S=0) | 1 | 8 |
| chi_c0 | a1 | rho | (L=1, S=1) | (L=0, S=1), (L=2, S=1) | (L=1, S=0) | 2 | 16 |
| chi_c0 | a1 | f2 | (L=1, S=1) | (L=1, S=2), (L=3, S=2) | (L=2, S=0) | 2 | 16 |
| chi_c0 | a2 | f0 | forbidden | forbidden | (L=0, S=0) | 0 | 0 |
| chi_c0 | a2 | rho | forbidden | (L=2, S=1) | (L=1, S=0) | 0 | 0 |
| chi_c0 | a2 | f2 | forbidden | (L=1, S=2), (L=3, S=2) | (L=2, S=0) | 0 | 0 |
| chi_c0 | pi0 | f0 | (L=0, S=0) | (L=0, S=0) | (L=0, S=0) | 1 | 8 |
| chi_c0 | pi0 | rho | (L=0, S=0) | (L=1, S=1) | (L=1, S=0) | 1 | 8 |
| chi_c0 | pi0 | f2 | (L=0, S=0) | (L=2, S=2) | (L=2, S=0) | 1 | 8 |
| chi_c0 | pi1 | f0 | forbidden | forbidden | (L=0, S=0) | 0 | 0 |
| chi_c0 | pi1 | rho | forbidden | (L=1, S=1) | (L=1, S=0) | 0 | 0 |
| chi_c0 | pi1 | f2 | forbidden | (L=2, S=2) | (L=2, S=0) | 0 | 0 |
| chi_c0 | pi2 | f0 | (L=2, S=2) | (L=2, S=0) | (L=0, S=0) | 1 | 8 |
| chi_c0 | pi2 | rho | (L=2, S=2) | (L=1, S=1), (L=3, S=1) | (L=1, S=0) | 2 | 16 |
| chi_c0 | pi2 | f2 | (L=2, S=2) | (L=0, S=2), (L=2, S=2), (L=4, S=2) | (L=2, S=0) | 3 | 24 |
| chi_c2 | a1 | f0 | (L=1, S=1), (L=3, S=1) | (L=1, S=0) | (L=0, S=0) | 2 | 16 |
| chi_c2 | a1 | rho | (L=1, S=1), (L=3, S=1) | (L=0, S=1), (L=2, S=1) | (L=1, S=0) | 4 | 32 |
| chi_c2 | a1 | f2 | (L=1, S=1), (L=3, S=1) | (L=1, S=2), (L=3, S=2) | (L=2, S=0) | 4 | 32 |
| chi_c2 | a2 | f0 | (L=1, S=2), (L=3, S=2) | forbidden | (L=0, S=0) | 0 | 0 |
| chi_c2 | a2 | rho | (L=1, S=2), (L=3, S=2) | (L=2, S=1) | (L=1, S=0) | 2 | 16 |
| chi_c2 | a2 | f2 | (L=1, S=2), (L=3, S=2) | (L=1, S=2), (L=3, S=2) | (L=2, S=0) | 4 | 32 |
| chi_c2 | pi0 | f0 | (L=2, S=0) | (L=0, S=0) | (L=0, S=0) | 1 | 8 |
| chi_c2 | pi0 | rho | (L=2, S=0) | (L=1, S=1) | (L=1, S=0) | 1 | 8 |
| chi_c2 | pi0 | f2 | (L=2, S=0) | (L=2, S=2) | (L=2, S=0) | 1 | 8 |
| chi_c2 | pi1 | f0 | (L=2, S=1) | forbidden | (L=0, S=0) | 0 | 0 |
| chi_c2 | pi1 | rho | (L=2, S=1) | (L=1, S=1) | (L=1, S=0) | 1 | 8 |
| chi_c2 | pi1 | f2 | (L=2, S=1) | (L=2, S=2) | (L=2, S=0) | 1 | 8 |
| chi_c2 | pi2 | f0 | (L=0, S=2), (L=2, S=2), (L=4, S=2) | (L=2, S=0) | (L=0, S=0) | 3 | 24 |
| chi_c2 | pi2 | rho | (L=0, S=2), (L=2, S=2), (L=4, S=2) | (L=1, S=1), (L=3, S=1) | (L=1, S=0) | 6 | 48 |
| chi_c2 | pi2 | f2 | (L=0, S=2), (L=2, S=2), (L=4, S=2) | (L=0, S=2), (L=2, S=2), (L=4, S=2) | (L=2, S=0) | 9 | 72 |

| parent | r2a | r2b | parent -\> R2a R2b | R2a -\> pi pi | R2b -\> pi pi | n_ls | explicit_chains |
|----|----|----|----|----|----|----|----|
| chi_c0 | f0 | f0 | (L=0, S=0) | (L=0, S=0) | (L=0, S=0) | 1 | 2 |
| chi_c0 | f0 | rho | (L=1, S=1) | (L=0, S=0) | (L=1, S=0) | 1 | 2 |
| chi_c0 | f0 | f2 | (L=2, S=2) | (L=0, S=0) | (L=2, S=0) | 1 | 2 |
| chi_c0 | rho | f0 | (L=1, S=1) | (L=1, S=0) | (L=0, S=0) | 1 | 2 |
| chi_c0 | rho | rho | (L=0, S=0), (L=2, S=2) | (L=1, S=0) | (L=1, S=0) | 2 | 4 |
| chi_c0 | rho | f2 | (L=1, S=1), (L=3, S=3) | (L=1, S=0) | (L=2, S=0) | 2 | 4 |
| chi_c0 | f2 | f0 | (L=2, S=2) | (L=2, S=0) | (L=0, S=0) | 1 | 2 |
| chi_c0 | f2 | rho | (L=1, S=1), (L=3, S=3) | (L=2, S=0) | (L=1, S=0) | 2 | 4 |
| chi_c0 | f2 | f2 | (L=0, S=0), (L=2, S=2), (L=4, S=4) | (L=2, S=0) | (L=2, S=0) | 3 | 6 |
| chi_c2 | f0 | f0 | (L=2, S=0) | (L=0, S=0) | (L=0, S=0) | 1 | 2 |
| chi_c2 | f0 | rho | (L=1, S=1), (L=3, S=1) | (L=0, S=0) | (L=1, S=0) | 2 | 4 |
| chi_c2 | f0 | f2 | (L=0, S=2), (L=2, S=2), (L=4, S=2) | (L=0, S=0) | (L=2, S=0) | 3 | 6 |
| chi_c2 | rho | f0 | (L=1, S=1), (L=3, S=1) | (L=1, S=0) | (L=0, S=0) | 2 | 4 |
| chi_c2 | rho | rho | (L=0, S=2), (L=2, S=0), (L=2, S=1), (L=2, S=2), (L=4, S=2) | (L=1, S=0) | (L=1, S=0) | 5 | 10 |
| chi_c2 | rho | f2 | (L=1, S=1), (L=1, S=2), (L=1, S=3), (L=3, S=1), (L=3, S=2), (L=3, S=3), (L=5, S=3) | (L=1, S=0) | (L=2, S=0) | 7 | 14 |
| chi_c2 | f2 | f0 | (L=0, S=2), (L=2, S=2), (L=4, S=2) | (L=2, S=0) | (L=0, S=0) | 3 | 6 |
| chi_c2 | f2 | rho | (L=1, S=1), (L=1, S=2), (L=1, S=3), (L=3, S=1), (L=3, S=2), (L=3, S=3), (L=5, S=3) | (L=2, S=0) | (L=1, S=0) | 7 | 14 |
| chi_c2 | f2 | f2 | (L=0, S=2), (L=2, S=0), (L=2, S=1), (L=2, S=2), (L=2, S=3), (L=2, S=4), (L=4, S=2), (L=4, S=3), (L=4, S=4), (L=6, S=4) | (L=2, S=0) | (L=2, S=0) | 10 | 20 |

## Dynamics

The dynamic model below is deliberately small. It uses the `chi_c2`
parent and the `a1 rho` cascade family, expanded over all eight
charge/Bose topology copies and all allowed LS products. The lineshapes
are toy Breit-Wigners; the numbers only make the benchmark finite and
reproducible.

| parent | stage                       | bracket            | r3  | r2  | n_chains |
|--------|-----------------------------|--------------------|-----|-----|----------|
| chi_c2 | positive: identity          | `(((1, 2), 3), 4)` | a1  | rho | 4        |
| chi_c2 | positive: swap pi+          | `(((3, 2), 1), 4)` | a1  | rho | 4        |
| chi_c2 | positive: swap pi-          | `(((1, 4), 3), 2)` | a1  | rho | 4        |
| chi_c2 | positive: swap pi+ and pi-  | `(((3, 4), 1), 2)` | a1  | rho | 4        |
| chi_c2 | conjugate: identity         | `(((2, 1), 4), 3)` | a1  | rho | 4        |
| chi_c2 | conjugate: swap pi+         | `(((2, 3), 4), 1)` | a1  | rho | 4        |
| chi_c2 | conjugate: swap pi-         | `(((4, 1), 2), 3)` | a1  | rho | 4        |
| chi_c2 | conjugate: swap pi+ and pi- | `(((4, 3), 2), 1)` | a1  | rho | 4        |

## Phase Space and Timing

The timing cells use only 100 flat phase-space events. The goal is a
local sanity check: how much time is spent turning four-vectors into
topology-local kinematics, and how much time is spent evaluating the toy
coherent amplitude.

| step                        | events | total_ms | ms_per_event |
|-----------------------------|--------|----------|--------------|
| KinematicPoint construction | 100    | 1228.34  | 12.2834      |
| coherent amplitude          | 100    | 1598.03  | 15.9803      |
| sum(abs2, amplitude)        | 100    | 29.2105  | 0.292105     |

The benchmark model is intentionally not the full physics model. It
exercises the same kinematic routing and amplitude machinery on a
compact, auditable subset. The full catalogue above tells which
additional resonance and LS families should be added when moving from
this toy benchmark to the analysis model.
