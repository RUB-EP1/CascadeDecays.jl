# Quarto is run locally. This file only adapts committed GFM for Documenter.
# Paths in `qmd` are relative to docs/; `page` paths are relative to the build source.
const RENDERED_PAGES = [
    (
        qmd = "kinematic-task.qmd", page = "kinematic-task.md",
        title = "Routing four-vectors", id = "kinematic_tasks",
    ),
    (
        qmd = "integration_4body_b2ddKpi.qmd", page = "tutorial.md",
        title = "Using a decay chain", id = "tutorial",
    ),
    (
        qmd = "lb2lc3pi-model.qmd", page = "lb2lc3pi-model.md",
        title = "Building a full model for a decay", id = "lb2lc3pi_model",
    ),
    (
        qmd = "cascade-vs-dpd.qmd", page = "cascade-vs-dpd.md",
        title = "Cross-checking with ThreeBodyDecays", id = "cascade_vs_dpd",
    ),
    (
        qmd = "particle-two-phase.qmd", page = "particle-two-phase.md",
        title = "Why do we need the particle-2 phase?", id = "particle_two_phase",
    ),
    (
        qmd = "four-pion-model.qmd", page = "four-pion-model.md",
        title = "Four-pion model-building catalogue", id = "four_pion_model",
    ),
    (
        qmd = "isospin-kaon-conventions.qmd", page = "isospin-kaon-conventions.md",
        title = "Isospin and kaon charge-conjugation conventions", id = "isospin_kaon_conventions",
    ),
    (
        qmd = "compass-basis.qmd", page = "compass-basis.md",
        title = "The COMPASS basis for three-pion amplitudes", id = "compass_basis",
    ),
    (
        qmd = "three-pion-cascade-vs-dpd.qmd", page = "three-pion-cascade-vs-dpd.md",
        title = "Three pions: pulling out the (12)3 orientation", id = "three_pion_cascade_vs_dpd",
    ),
    (
        qmd = "../examples/pp2ppKK-model.qmd", page = "examples/pp2ppKK-model.md",
        title = "A pp to pp K+ K- amplitude model", id = "pp2ppKK_model",
    ),
]

function documenter_math_blocks(body::AbstractString)
    lines = String.(split(body, '\n'))
    fence = ""
    in_math = false
    for (i, line) in pairs(lines)
        marker = match(r"^ {0,3}(`{3,}|~{3,})(.*)$", line)
        if !isempty(fence)
            if marker !== nothing && startswith(marker[1], fence) && isempty(strip(marker[2]))
                fence = ""
            end
        elseif in_math
            if strip(line) == raw"$$"
                lines[i] = "```"
                in_math = false
            end
        elseif marker !== nothing
            fence = marker[1]
        elseif strip(line) == raw"$$"
            lines[i] = "```math"
            in_math = true
        end
    end
    in_math && error("Unclosed display-math block in rendered Markdown")
    return join(lines, '\n')
end

function documenter_quarto_page(body::AbstractString; edit_url, title, id)
    # Julia Markdown needs fenced math for multiline GFM display equations.
    body = documenter_math_blocks(body)
    heading = "# [$title](@id $id)"
    body = replace(body, r"^# [^\r\n]+\r?\n" => heading * "\n"; count = 1)
    # Pandoc URL-encodes the space in named Documenter cross-references.
    body = replace(body, "(@ref%20" => "(@ref ")
    body = replace(body, r"<img\s+src=\"([^\"]+)\"\s+id=\"([^\"]+)\"\s*/>" => s"![](\1)")
    meta = "```@meta\nCurrentModule = CascadeDecays\nEditURL = \"$edit_url\"\n```\n\n"
    return meta * body
end

function copy_rendered_page!(spec, docs, destination_root)
    qmd = normpath(joinpath(docs, spec.qmd))
    gfm = joinpath(dirname(qmd), "generated", replace(basename(qmd), r"\.qmd$" => ".md"))
    command_path = replace(relpath(qmd, dirname(docs)), '\\' => '/')
    hint = "Run locally: quarto render $command_path --to gfm; commit the Markdown and figures."
    isfile(gfm) || error("Missing rendered Markdown: $gfm. $hint")
    body = read(gfm, String)
    assets = splitext(basename(gfm))[1] * "_files"
    source_assets = joinpath(dirname(gfm), assets)
    occursin(assets * "/", body) && !isdir(source_assets) &&
        error("Missing rendered figures: $source_assets. $hint")

    destination = joinpath(destination_root, spec.page)
    mkpath(dirname(destination))
    # EditURL is interpreted relative to the temporary Documenter source page.
    edit_url = replace(relpath(qmd, dirname(destination)), '\\' => '/')
    write(destination, documenter_quarto_page(body; edit_url, spec.title, spec.id))
    if occursin(assets * "/", body)
        cp(source_assets, joinpath(dirname(destination), assets); force = true)
    end
    return destination
end

function prepare_documenter_source!(;
        docs = @__DIR__, destination = joinpath(docs, ".documenter-src"),
    )
    rm(destination; recursive = true, force = true)
    mkpath(destination)
    # Ignore old generated copies left by earlier builds in docs/src/.
    generated_names = Set(basename(spec.page) for spec in RENDERED_PAGES)
    union!(generated_names, (splitext(basename(spec.qmd))[1] * "_files" for spec in RENDERED_PAGES))
    for entry in readdir(joinpath(docs, "src"))
        entry in generated_names && continue
        source = joinpath(docs, "src", entry)
        isdir(source) || endswith(entry, ".md") || continue
        cp(source, joinpath(destination, entry); force = true)
    end
    foreach(spec -> copy_rendered_page!(spec, docs, destination), RENDERED_PAGES)
    return destination
end
