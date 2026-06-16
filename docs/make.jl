using CascadeDecays
using Documenter

const DOCS = @__DIR__
const TUTORIAL_QMD = joinpath(DOCS, "integration_4body_b2ddKpi.qmd")
const TUTORIAL_GFM = joinpath(DOCS, "integration_4body_b2ddKpi.md")
const TUTORIAL = joinpath(DOCS, "src", "tutorial.md")
const LB2LC3PI_QMD = joinpath(DOCS, "lb2lc3pi-model.qmd")
const LB2LC3PI_GFM = joinpath(DOCS, "lb2lc3pi-model.md")
const LB2LC3PI_PAGE = joinpath(DOCS, "src", "lb2lc3pi-model.md")

function render_quarto_gfm!(qmd_path::AbstractString, gfm_path::AbstractString)
    cd(DOCS) do
        run(`quarto render $(basename(qmd_path)) --to gfm`)
    end
    isfile(gfm_path) || error("expected Quarto output at $(gfm_path)")
end

function copy_quarto_assets!(gfm_path::AbstractString, page_path::AbstractString)
    asset_dir = splitext(basename(gfm_path))[1] * "_files"
    source = joinpath(dirname(gfm_path), asset_dir)
    destination = joinpath(dirname(page_path), asset_dir)
    isdir(source) || return nothing
    isdir(destination) && rm(destination; recursive = true)
    cp(source, destination)
    return nothing
end

function documenter_quarto_page(gfm_path::AbstractString; edit_url::AbstractString, title = nothing)
    body = read(gfm_path, String)
    if title !== nothing
        body = replace(body, r"^# .+\n\n" => title * "\n\n"; count = 1)
    end
    body = replace(body, r"<img\s+src=\"([^\"]+)\"\s+id=\"([^\"]+)\"\s*/>" => s"![](\1)")
    meta = "```@meta\nCurrentModule = CascadeDecays\nEditURL = \"$(edit_url)\"\n```\n\n"
    return meta * body
end

DocMeta.setdocmeta!(
    CascadeDecays,
    :DocTestSetup,
    :(using CascadeDecays);
    recursive = true,
)

render_quarto_gfm!(TUTORIAL_QMD, TUTORIAL_GFM)
write(TUTORIAL, documenter_quarto_page(TUTORIAL_GFM; edit_url = "../integration_4body_b2ddKpi.qmd"))

render_quarto_gfm!(LB2LC3PI_QMD, LB2LC3PI_GFM)
copy_quarto_assets!(LB2LC3PI_GFM, LB2LC3PI_PAGE)
write(
    LB2LC3PI_PAGE,
    documenter_quarto_page(
        LB2LC3PI_GFM;
        edit_url = "../lb2lc3pi-model.qmd",
        title = "# [``\\Lambda_b^0 \\to \\Lambda_c^+ \\pi^+ \\pi^- \\pi^-`` model setup](@id lb2lc3pi_model)",
    ),
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
        "Internal notation" => "notation.md",
        "Kinematic tasks" => "kinematic-task.md",
        "Tutorial" => "tutorial.md",
        "Lb to Lc 3pi model" => "lb2lc3pi-model.md",
        "CascadeDecays vs DPD" => "cascade-vs-dpd.md",
        "API reference" => "api-reference.md",
    ],
)

deploydocs(; repo = "github.com/RUB-EP1/CascadeDecays.jl")
