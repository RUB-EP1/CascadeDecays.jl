# CascadeDecays examples

This directory keeps executable Quarto notebooks together with their committed,
pre-rendered Markdown and plot assets. Documenter publishes the rendered files;
it does not execute these notebooks or regenerate their Monte Carlo samples in
CI.

## Layout

Each example consists of:

- `name.qmd`: the executable source notebook;
- `name.md`: the committed GitHub-flavored Markdown render;
- `name_files/`: committed figures and other render assets.

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
git diff -- examples/pp2ppKK-model.qmd examples/pp2ppKK-model.md
git status --short examples/pp2ppKK-model_files
```

To verify that Documenter can consume the static render without re-running the
sample, run the documentation build normally:

```sh
julia --project=docs docs/make.jl
```

The documentation build checks that the committed Markdown and asset directory
exist, copies them into Documenter's temporary source tree, and fails with a
regeneration hint if either is missing.
