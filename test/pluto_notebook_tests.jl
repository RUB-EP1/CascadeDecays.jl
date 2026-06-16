# Pluto notebook structure tests (parse-only; no ServerSession).
# Requires Pluto in the active environment: ] add Pluto

using Test

@testset "Pluto notebook structure" begin
    using Pluto

    notebook = joinpath(@__DIR__, "..", "docs", "Lc2pKpi_L_Kstar.jl")
    @test isfile(notebook)
    @test readline(notebook) == "### A Pluto.jl notebook ###"
    @test Pluto.is_pluto_notebook(notebook)

    nb = Pluto.load_notebook_nobackup(notebook)
    @test length(nb.cells) >= 1
    @test length(nb.cell_order) == length(nb.cells)
    @test all(haskey(nb.cells_dict, id) for id in nb.cell_order)
end
