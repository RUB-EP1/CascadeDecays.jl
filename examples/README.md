# CascadeDecays examples

This directory keeps executable Quarto sources and stores their committed
Markdown and plot assets in `generated/`. Documenter publishes the rendered files;
it does not execute these notebooks or regenerate their Monte Carlo samples in
CI.

## Layout

Each example consists of:

- `name.qmd`: the executable source notebook;
- `generated/name.md`: the committed GitHub-flavored Markdown render;
- `generated/name_files/`: committed figures and other render assets.

Keep all three in sync when changing an example.

## Regenerating an example

The notebooks use the Julia environment in this directory. From the repository
root, instantiate it once:

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

Then render the requested notebook from the repository root:

```sh
quarto render examples/pp2ppKK-model.qmd --to gfm
```

Review and commit the source, rendered Markdown, and generated asset directory:

```sh
git diff -- examples/pp2ppKK-model.qmd examples/generated/pp2ppKK-model.md
git status --short examples/generated/pp2ppKK-model_files
```

To verify that Documenter can consume the static render without re-running the
sample, run the documentation build normally:

```sh
julia --project=docs docs/make.jl
```

The documentation build checks that the committed Markdown and asset directory
exist, copies them into Documenter's temporary source tree, and fails with a
regeneration hint if either is missing.

## Four-body sampling in a three-pion mass window

`validate-fourbody-mass-window.jl` is a standalone numerical check (no notebook
or generated assets). From the repository root, using the examples environment:

```sh
julia -O1 --project=examples examples/validate-fourbody-mass-window.jl
```

The `-O1` flag avoids an LLVM optimizer crash observed with Julia 1.11.5 on
Apple Silicon. The examples environment pins RamboOnDiet `v0.5.0`, the first
release with the massive-weight normalization correction. Earlier versions
give the right shape but the wrong absolute normalization, which this script
tests.

The test generates 100,000 events for a 190 GeV/c pion beam on a stationary
proton, ordering outputs as `[recoil, pi-, pi+, pi-]`. Restricting the first
cube coordinate places every three-pion mass in `[1.02,1.7]` GeV. Event weights
include the width of that coordinate interval. An optional first argument
changes the sample size (minimum 10,000).

It checks mass-window endpoints, conservation, mass shells, the raw sampling
distribution, and the weighted spectrum against independent two-body
phase-space factorization with Gauss-Legendre quadrature. The phase-space
convention is `delta^4(P-sum(p)) product(d^3p/(2E))`, without `2pi` factors.
The quadrature is repeated at two orders to check convergence.

With seed 20261021, all 20 sampling and kinematic assertions plus two
serialization round-trip checks pass. The restricted volume is
`5.1253743 +/- 0.00422`, versus quadrature `5.1311734` (1.37 standard errors).
The effective weighted sample is 93.65%; the largest of ten mass-bin
residuals is 1.57 standard errors. These are weighted phase-space events,
not an unweighted massive phase-space sample.

The script saves `examples/data/compass-mass-sample.jls`, a Julia-serialized
named tuple containing `metadata` and column arrays in `columns`:

| Column | Definition |
|---|---|
| `m3pi`, `m12`, `m23` | Masses in GeV; labels 1=pi-, 2=pi+, 3=pi- |
| `theta12` | TY polar angle of particle 1 in pair (12), radians |
| `theta23` | TY polar angle of particle 2 in pair (23), radians |
| `zeta0_13` | Signed negative parent rotation `wr(1,3,0)` from (12) reference to cyclic (23), radians |
| `weight` | Restricted four-body integration weight, including cube-window width |

No four-vectors or orientation angles are stored. The three stored polar
angles are reconstructed from the three masses and checked independently
against Cartesian boosts for every event (maximum discrepancy
`6.36e-11` radians). The regenerated sample uses schema version 2: `m13`, `theta31`, and
`zeta0_23` have been replaced by `m23`, `theta23`, and signed `zeta0_13`.
The metadata records masses, beam momentum, seed,
window, units, ordering, measure, and weight convention. This is a native
Julia binary format; regenerate it using the script when changing Julia versions.

```julia
using Serialization, Statistics
sample = deserialize("examples/data/compass-mass-sample.jls")
x = sample.columns
# For a vector f containing one integrand value per event:
# integral = mean(x.weight .* f)
```

The weight estimates the restricted **four-body** measure after integrating
unused orientations; it does not convert the sample to `dm3pi * dPhi3`.
The generated binary sample is ignored by Git.
