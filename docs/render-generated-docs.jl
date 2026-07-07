const DOCS = @__DIR__

const GENERATED_PAGES = [
    (
        qmd = "four-pion-model.qmd",
        gfm = "four-pion-model.md",
        dest = joinpath("generated", "four-pion-model.md"),
        title = "Four-pion model-building catalogue",
        id = "four_pion_model",
    ),
    (
        qmd = "isospin-kaon-conventions.qmd",
        gfm = "isospin-kaon-conventions.md",
        dest = joinpath("generated", "isospin-kaon-conventions.md"),
        title = "Isospin and kaon charge-conjugation conventions",
        id = "isospin_kaon_conventions",
    ),
]

function documenter_page(body::AbstractString; title, id, qmd)
    body = replace(
        body,
        r"^# .+\r?\n\r?\n" => "# [$title](@id $id)\n\n";
        count = 1,
    )
    meta = "```@meta\nCurrentModule = CascadeDecays\nEditURL = \"../$qmd\"\n```\n\n"
    return meta * body
end

function render_generated_page!(page)
    cd(DOCS) do
        run(`quarto render $(page.qmd) --to gfm`)
    end
    gfm_path = joinpath(DOCS, page.gfm)
    isfile(gfm_path) || error("expected Quarto output at $(gfm_path)")
    dest_path = joinpath(DOCS, page.dest)
    write(
        dest_path,
        documenter_page(read(gfm_path, String); title = page.title, id = page.id, qmd = page.qmd),
    )
    rm(gfm_path; force = true)
end

for page in GENERATED_PAGES
    render_generated_page!(page)
end
