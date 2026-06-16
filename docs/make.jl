using CascadeDecays
using Documenter

const DOCS = @__DIR__
const QMD = joinpath(DOCS, "integration_4body_b2ddKpi.qmd")
const GFM = joinpath(DOCS, "integration_4body_b2ddKpi.md")
const TUTORIAL = joinpath(DOCS, "src", "tutorial.md")

function render_integration_tutorial!()
    cd(DOCS) do
        run(`quarto render $(basename(QMD)) --to gfm`)
    end
    isfile(GFM) || error("expected Quarto output at $(GFM)")
end

function documenter_tutorial_page(gfm_path::AbstractString)
    body = read(gfm_path, String)
    meta = "```@meta\nCurrentModule = CascadeDecays\nEditURL = \"../integration_4body_b2ddKpi.qmd\"\n```\n\n"
    return meta * body
end

DocMeta.setdocmeta!(
    CascadeDecays,
    :DocTestSetup,
    :(using CascadeDecays);
    recursive = true,
)

render_integration_tutorial!()
write(TUTORIAL, documenter_tutorial_page(GFM))

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
        "Internal notation" => "notation.md",
        "Kinematic tasks" => "kinematic-task.md",
        "Tutorial" => "tutorial.md",
        "CascadeDecays vs DPD" => "cascade-vs-dpd.md",
        "API reference" => "api-reference.md",
    ],
)

deploydocs(; repo = "github.com/RUB-EP1/CascadeDecays.jl")
