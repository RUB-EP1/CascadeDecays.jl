#!/usr/bin/env julia
# Validate Documenter cross-references after `docs/make.jl`.

const DOCS = @__DIR__
const BUILD = joinpath(DOCS, "build")

function page_ids()
    return Dict(
        "home" => "index.html",
        "notation" => "notation/index.html",
        "kinematic_tasks" => "kinematic-task/index.html",
        "tutorial" => "tutorial/index.html",
        "lb2lc3pi_model" => "lb2lc3pi-model/index.html",
        "cascade_vs_dpd" => "cascade-vs-dpd/index.html",
        "api_reference" => "api-reference/index.html",
    )
end

function collect_hrefs(html::AbstractString)
    hrefs = String[]
    for m in eachmatch(r"href=\"([^\"]+)\"", html)
        push!(hrefs, m.captures[1])
    end
    return hrefs
end

function validate_page_ids!(ids::Dict{String,String})
    for (id, relpath) in ids
        path = joinpath(BUILD, relpath)
        isfile(path) || error("missing built page for @id $id: $path")
        html = read(path, String)
        occursin("#$id", html) || occursin("id=\"$id\"", html) ||
            error("built page for @id $id does not contain anchor #$id: $relpath")
    end
end

function validate_doc_map_links!()
    index = read(joinpath(BUILD, "index.html"), String)
    expected = [
        "notation/#notation",
        "kinematic-task/#kinematic_tasks",
        "tutorial/#tutorial",
        "lb2lc3pi-model/#lb2lc3pi_model",
        "cascade-vs-dpd/#cascade_vs_dpd",
        "api-reference/#api_reference",
    ]
    for target in expected
        occursin(target, index) || error("documentation map missing link target: $target")
    end
end

function validate_no_broken_markers!()
    for (root, _, files) in walkdir(BUILD)
        for file in files
            endswith(file, ".html") || continue
            html = read(joinpath(root, file), String)
            if occursin("Reference not found", html) || occursin("@ref", html)
                error("possible unresolved @ref in $(joinpath(root, file))")
            end
        end
    end
end

function main()
    isdir(BUILD) || error("missing docs/build; run `julia --project=docs docs/make.jl` first")
    validate_page_ids!(page_ids())
    validate_doc_map_links!()
    validate_no_broken_markers!()
    println("Cross-reference validation passed.")
end

main()
