const DOCS = @__DIR__
const SOURCE_DIR = joinpath(DOCS, "isospin_kaon_note_updated")
const QMD = joinpath(SOURCE_DIR, "isospin-kaon-conventions.qmd")
const GFM = joinpath(SOURCE_DIR, "isospin-kaon-conventions.md")
const DEST = joinpath(DOCS, "src", "isospin-kaon-conventions.md")
const TITLE = "Isospin and kaon charge-conjugation conventions"

function strip_frontmatter(text::AbstractString)
    lines = split(text, '\n'; keepempty=true)
    length(lines) >= 2 && lines[1] == "---" || return text
    close = findnext(==("---"), lines, 2)
    close === nothing && return text
    return join(lines[(close + 1):end], '\n')
end

function source_artifact(body::AbstractString)
    return "# $TITLE\n\n" * body
end

function documenter_page(body::AbstractString)
    meta = "```@meta\nCurrentModule = CascadeDecays\nEditURL = \"../isospin_kaon_note_updated/isospin-kaon-conventions.qmd\"\n```\n\n"
    header = "# [$TITLE](@id isospin_kaon_conventions)\n\n"
    return meta * header * body
end

body = strip_frontmatter(read(QMD, String))
write(GFM, source_artifact(body))
write(DEST, documenter_page(body))
