using CascadeDecays
using Documenter

const DOCS = @__DIR__

# Quarto notebooks rendered into Documenter pages under src/.
const QUARTO_PAGES = [
    (
        qmd = "integration_4body_b2ddKpi.qmd",
        page = "tutorial.md",
        edit_url = "../integration_4body_b2ddKpi.qmd",
    ),
    (
        qmd = "lb2lc3pi-model.qmd",
        page = "lb2lc3pi-model.md",
        edit_url = "../lb2lc3pi-model.qmd",
        title = "# [``\\Lambda_b^0 \\to \\Lambda_c^+ \\pi^+ \\pi^- \\pi^-`` model setup](@id lb2lc3pi_model)",
        copy_assets = true,
    ),
]

gfm_path(qmd::AbstractString) = joinpath(DOCS, replace(qmd, r"\.qmd$" => ".md"))
page_path(page::AbstractString) = joinpath(DOCS, "src", page)

function render_quarto_gfm!(qmd::AbstractString)
    cd(DOCS) do
        run(`quarto render $qmd --to gfm`)
    end
    path = gfm_path(qmd)
    isfile(path) || error("expected Quarto output at $(path)")
    return path
end

function copy_quarto_assets!(gfm_path::AbstractString, page_path::AbstractString)
    asset_dir = splitext(basename(gfm_path))[1] * "_files"
    source = joinpath(dirname(gfm_path), asset_dir)
    destination = joinpath(dirname(page_path), asset_dir)
    isdir(source) || return
    isdir(destination) && rm(destination; recursive = true)
    cp(source, destination)
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

function build_quarto_page!(spec)
    gfm = render_quarto_gfm!(spec.qmd)
    dest = page_path(spec.page)
    get(spec, :copy_assets, false) && copy_quarto_assets!(gfm, dest)
    write(
        dest,
        documenter_quarto_page(gfm; edit_url = spec.edit_url, title = get(spec, :title, nothing)),
    )
end

function build_quarto_pages!()
    foreach(build_quarto_page!, QUARTO_PAGES)
end

DocMeta.setdocmeta!(
    CascadeDecays,
    :DocTestSetup,
    :(using CascadeDecays);
    recursive = true,
)

build_quarto_pages!()

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
