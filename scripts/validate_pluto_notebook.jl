#!/usr/bin/env julia
# Validate Pluto notebook file structure (and optionally run as a script).
#
# Structure check without Pluto (header + cell-order footer only):
#   julia scripts/validate_pluto_notebook.jl --header-only docs/Lc2pKpi_L_Kstar.jl
#
# Full parse via Pluto.load_notebook_nobackup (no cell execution):
#   julia --project=. scripts/validate_pluto_notebook.jl
#   julia --project=. scripts/validate_pluto_notebook.jl docs/Lc2pKpi_L_Kstar.jl
#
# Also execute top-to-bottom as a plain Julia script (mock @bind defaults):
#   julia --project=. scripts/validate_pluto_notebook.jl --run docs/Lc2pKpi_L_Kstar.jl

const NOTEBOOK_HEADER = "### A Pluto.jl notebook ###"
const DEFAULT_NOTEBOOK = joinpath(@__DIR__, "..", "docs", "Lc2pKpi_L_Kstar.jl")

if !any(==("--header-only"), ARGS)
    using Pluto
end

function parse_args(ARGS)
    run_script = false
    header_only = false
    paths = String[]
    for arg in ARGS
        if arg == "--run"
            run_script = true
        elseif arg == "--header-only"
            header_only = true
        elseif startswith(arg, "-")
            error("Unknown option: $arg")
        else
            push!(paths, abspath(arg))
        end
    end
    isempty(paths) && push!(paths, abspath(DEFAULT_NOTEBOOK))
    return (; run_script, header_only, paths)
end

"""Cheap check: Pluto header, .jl extension, cell-order block present."""
function validate_header_only(path::AbstractString)
    path = abspath(path)
    isfile(path) || error("File not found: $path")
    readline(path) == NOTEBOOK_HEADER ||
        error("Invalid header (expected Pluto notebook marker): $path")
    endswith(path, ".jl") || error("Expected a .jl notebook path: $path")
    text = read(path, String)
    occursin("# ╔═╡ Cell order:", text) ||
        error("Missing cell-order footer: $path")
    n_delims = count(line -> startswith(line, "# ╔═╡ "), split(text, '\n')) - 1
    n_delims >= 1 || error("No cells found: $path")
    return (; path, n_cells = n_delims)
end

# Parse-only validation using Pluto.load_notebook_nobackup (no cell execution).
function validate_structure(path::AbstractString)
    path = abspath(path)
    isfile(path) || error("File not found: $path")
    first_line = readline(path)
    if first_line != "### A Pluto.jl notebook ###"
        error("Invalid header (expected Pluto notebook marker): $path")
    end
    if !Pluto.is_pluto_notebook(path)
        error("Pluto.is_pluto_notebook returned false: $path")
    end
    nb = Pluto.load_notebook_nobackup(path)
    n_cells = length(nb.cells)
    n_cells >= 1 || error("Notebook has no cells: $path")
    length(nb.cell_order) == n_cells ||
        error("cell_order length ($(length(nb.cell_order))) != cells ($n_cells): $path")
    for id in nb.cell_order
        haskey(nb.cells_dict, id) || error("cell_order UUID missing from cells_dict: $id")
    end
    return (; path, n_cells, cell_order = nb.cell_order)
end

function validate_run(path::AbstractString)
    path = abspath(path)
    project = joinpath(@__DIR__, "..")
    proc = run(ignorestatus(setenv(`$(Base.julia_cmd()) --project=$project $path`; dir = dirname(path))))
    proc.exitcode == 0 || error("Script run failed (exit $(proc.exitcode)): $path")
    return path
end

function main(ARGS)
    opts = parse_args(ARGS)
    for path in opts.paths
        if opts.header_only
            info = validate_header_only(path)
            println("OK header     ", info.path, "  (", info.n_cells, " cells)")
        else
            info = validate_structure(path)
            println("OK structure  ", info.path, "  (", info.n_cells, " cells)")
        end
        if opts.run_script
            validate_run(path)
            println("OK run        ", path)
        end
    end
end

main(ARGS)
