const DOCS = @__DIR__
const QMD = joinpath(DOCS, "isospin-kaon-conventions.qmd")
const DEST = joinpath(DOCS, "generated", "isospin-kaon-conventions.md")
const TITLE = "Isospin and kaon charge-conjugation conventions"

function strip_frontmatter(text::AbstractString)
    lines = split(text, '\n'; keepempty=true)
    length(lines) >= 2 && strip(lines[1]) == "---" || return text
    close = findnext(line -> strip(line) == "---", lines, 2)
    close === nothing && return text
    return join(lines[(close + 1):end], '\n')
end

function documenter_page(body::AbstractString)
    meta = "```@meta\nCurrentModule = CascadeDecays\nEditURL = \"../isospin-kaon-conventions.qmd\"\n```\n\n"
    header = "# [$TITLE](@id isospin_kaon_conventions)\n\n"
    return meta * header * body
end

body = strip_frontmatter(read(QMD, String))
write(DEST, documenter_page(body))
