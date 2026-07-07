const DOCS = @__DIR__
const QMD = "four-pion-model.qmd"
const GFM = joinpath(DOCS, "four-pion-model.md")
const DEST = joinpath(DOCS, "src", "four-pion-model.md")

function documenter_page(body::AbstractString)
    body = replace(
        body,
        r"^# .+\r?\n\r?\n" =>
            "# [Four-pion model-building catalogue](@id four_pion_model)\n\n";
        count = 1,
    )
    meta = "```@meta\nCurrentModule = CascadeDecays\nEditURL = \"../four-pion-model.qmd\"\n```\n\n"
    return meta * body
end

cd(DOCS) do
    run(`quarto render $QMD --to gfm`)
end

isfile(GFM) || error("expected Quarto output at $GFM")
write(DEST, documenter_page(read(GFM, String)))
