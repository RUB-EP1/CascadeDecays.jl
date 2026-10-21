# CascadeDecays.jl

[![Test workflow status](https://github.com/RUB-EP1/CascadeDecays.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/RUB-EP1/CascadeDecays.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Docs workflow status](https://github.com/RUB-EP1/CascadeDecays.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/RUB-EP1/CascadeDecays.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Dev documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://rub-ep1.github.io/CascadeDecays.jl/dev)

`CascadeDecays.jl` provides flat topology and routing utilities for sequential decay amplitudes. It connects static cascade definitions to event-by-event invariants and helicity-angle kinematics.

## Installation

With the checked-in `Manifest.toml` and Julia `1.11.5`, instantiate the project directly:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

If you work with a different Julia version than the one used for this manifest, resolve the environment explicitly:

```julia
using Pkg
Pkg.activate(".")
Pkg.resolve()
Pkg.instantiate()
```

## Development

- Run tests with `julia --project=. test/runtests.jl`.
- Render any documentation notebook locally with `quarto render docs/<name>.qmd --to gfm` (requires [Quarto](https://quarto.org/)). Commit the `.qmd`, sibling `.md`, and any `<name>_files/` figures together. See [docs/README.md](docs/README.md).
- Build the docs with `julia --project=docs -e 'using Pkg; Pkg.instantiate()'` followed by `julia --project=docs docs/make.jl`. Documenter consumes committed Markdown and resolves `@ref` links; it does not run Quarto or notebook computations, locally or in CI.
- When editing a tutorial, check its `execute.output` setting. For notebooks with `output: false`, add `#| output: true` to cells that should print tables or show return values; end other cells with `nothing` to avoid leaking the last expression.
- Handwritten pages live in `docs/src/`. Rendered notebooks are adapted for Documenter only in the ignored `docs/.documenter-src/` directory; do not edit those temporary copies.
- Before every major release, regenerate and review all documentation, tutorials, and examples using the [documentation maintenance checklist](docs/README.md#maintenance-before-major-releases) to catch problems and drift that the static CI build cannot detect.
