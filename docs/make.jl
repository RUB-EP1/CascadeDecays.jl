using CascadeDecays
using Documenter

include("prepare-pages.jl")

DocMeta.setdocmeta!(CascadeDecays, :DocTestSetup, :(using CascadeDecays); recursive = true)
prepare_documenter_source!()

makedocs(;
    modules = [CascadeDecays],
    authors = "Mikhail Mikhasenko and contributors",
    repo = "https://github.com/RUB-EP1/CascadeDecays.jl/blob/{commit}{path}#{line}",
    sitename = "CascadeDecays.jl",
    source = ".documenter-src",
    build = ".documenter-build",
    doctest = false,
    checkdocs = :none,
    format = Documenter.HTML(;
        canonical = "https://rub-ep1.github.io/CascadeDecays.jl",
        repolink = "https://github.com/RUB-EP1/CascadeDecays.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Topology and numbering" => "notation.md",
        "Amplitude computation" => "amplitude-computation.md",
        "Why do we need the particle-2 phase?" => "particle-two-phase.md",
        "Convention matching" => "convention-matching.md",
        "Routing four-vectors" => "kinematic-task.md",
        "Using a decay chain" => "tutorial.md",
        "Building a full model for a decay" => "lb2lc3pi-model.md",
        "Examples" => [
            "pp to pp K+ K- amplitude model" => "examples/pp2ppKK-model.md",
            "Four-pion model-building catalogue" => "four-pion-model.md",
            "Isospin and kaon charge-conjugation conventions" => "isospin-kaon-conventions.md",
        ],
        "Cross-checking with ThreeBodyDecays" => [
            "General relations and signs" => "cascade-vs-dpd.md",
            "COMPASS: three-pion meson amplitudes" => "compass-basis.md",
        ],
        hide("three-pion-cascade-vs-dpd.md"),
        "API reference" => "api-reference.md",
    ],
)

deploydocs(; repo = "github.com/RUB-EP1/CascadeDecays.jl", target = ".documenter-build")
