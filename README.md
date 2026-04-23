# CascadeDecays.jl

[![Test workflow status](https://github.com/mmikhasenko/CascadeDecays.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/CascadeDecays.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Docs workflow status](https://github.com/mmikhasenko/CascadeDecays.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/CascadeDecays.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Dev documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmikhasenko.github.io/CascadeDecays.jl/dev)

`CascadeDecays.jl` provides flat topology and routing utilities for sequential decay amplitudes. It connects static cascade definitions to event-by-event invariants and helicity-angle kinematics.

## Installation

With the checked-in `Manifest.toml` and Julia `1.10`, instantiate the project directly:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

If you work with a different Julia version than the one used for this manifest, resolve the environment explicitly and install the non-registered dependencies first:

```julia
using Pkg
Pkg.activate(".")
Pkg.develop(url = "https://github.com/mmikhasenko/FourVectors.jl")
Pkg.develop(url = "https://github.com/mmikhasenko/InstructionalDecayTrees.jl")
Pkg.instantiate()
```

## Development

- Run tests with `julia --project=. test/runtests.jl`.
- Build the docs with `julia --project=docs -e 'using Pkg; Pkg.instantiate()'` followed by `julia --project=docs docs/make.jl`.
- The minimal documentation includes the four-body integration tutorial in [docs/integration_4body_b2ddKpi.jl](/Users/mikhailmikhasenko/Documents/JuliaDev.CAT/CascadeDecays.jl/docs/integration_4body_b2ddKpi.jl).
