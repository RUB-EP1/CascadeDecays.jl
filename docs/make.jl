using CascadeDecays
using Documenter

const DOCS = @__DIR__
const DOCUMENTER_SOURCE = joinpath(DOCS, ".documenter-src")

# Quarto notebooks rendered into Documenter pages under src/.
const QUARTO_PAGES = [
    (
        qmd = "integration_4body_b2ddKpi.qmd",
        page = "tutorial.md",
        edit_url = "../integration_4body_b2ddKpi.qmd",
        title = "# [Using a decay chain](@id tutorial)",
    ),
    (
        qmd = "lb2lc3pi-model.qmd",
        page = "lb2lc3pi-model.md",
        edit_url = "../lb2lc3pi-model.qmd",
        title = "# [Building a full model for a decay](@id lb2lc3pi_model)",
        copy_assets = true,
    ),
    (
        qmd = "cascade-vs-dpd.qmd",
        page = "cascade-vs-dpd.md",
        edit_url = "../cascade-vs-dpd.qmd",
        title = "# [Cross-checking with ThreeBodyDecays](@id cascade_vs_dpd)",
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
    return cp(source, destination)
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
    return write(
        dest,
        documenter_quarto_page(gfm; edit_url = spec.edit_url, title = get(spec, :title, nothing)),
    )
end

function build_quarto_pages!()
    return foreach(build_quarto_page!, QUARTO_PAGES)
end

DocMeta.setdocmeta!(
    CascadeDecays,
    :DocTestSetup,
    :(using CascadeDecays);
    recursive = true,
)

build_quarto_pages!()

function prepare_documenter_source!()
    rm(DOCUMENTER_SOURCE; recursive = true, force = true)
    mkpath(DOCUMENTER_SOURCE)

    for entry in readdir(joinpath(DOCS, "src"))
        source = joinpath(DOCS, "src", entry)
        isdir(source) || splitext(entry)[2] == ".md" || continue
        cp(joinpath(DOCS, "src", entry), joinpath(DOCUMENTER_SOURCE, entry); force = true)
    end
    for entry in readdir(joinpath(DOCS, "generated"))
        cp(joinpath(DOCS, "generated", entry), joinpath(DOCUMENTER_SOURCE, entry); force = true)
    end
end

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
        "Routing four-vectors" => "kinematic-task.md",
        "Using a decay chain" => "tutorial.md",
        "Building a full model for a decay" => "lb2lc3pi-model.md",
        "Four-pion model-building catalogue" => "four-pion-model.md",
        "Isospin and kaon charge-conjugation conventions" => "isospin-kaon-conventions.md",
        "Cross-checking with ThreeBodyDecays" => "cascade-vs-dpd.md",
        "API reference" => "api-reference.md",
    ],
)

deploydocs(; repo = "github.com/RUB-EP1/CascadeDecays.jl", target = ".documenter-build")
