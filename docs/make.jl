using CascadeDecays
using Documenter
using Literate

DocMeta.setdocmeta!(
    CascadeDecays,
    :DocTestSetup,
    :(using CascadeDecays);
    recursive = true,
)

docs_src_dir = joinpath(@__DIR__, "src")
tutorial_src = joinpath(@__DIR__, "integration_4body_b2ddKpi.jl")

function strip_edit_url(content)
    replace(content, "EditURL = \"@__REPO_ROOT_URL__/\"" => "")
end

Literate.markdown(
    tutorial_src,
    docs_src_dir;
    name = "tutorial",
    documenter = true,
    credit = true,
    postprocess = strip_edit_url,
)

makedocs(;
    modules = [CascadeDecays],
    authors = "Mikhail Mikhasenko and contributors",
    repo = "https://github.com/RUB-EP1/CascadeDecays.jl/blob/{commit}{path}#{line}",
    sitename = "CascadeDecays.jl",
    doctest = false,
    checkdocs = :none,
    format = Documenter.HTML(;
        canonical = "https://rub-ep1.github.io/CascadeDecays.jl",
        repolink = "https://github.com/RUB-EP1/CascadeDecays.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
    ],
)

deploydocs(; repo = "github.com/RUB-EP1/CascadeDecays.jl")
