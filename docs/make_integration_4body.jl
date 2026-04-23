using Literate

const DOCS = @__DIR__
const INPUT = joinpath(DOCS, "integration_4body_b2ddKpi.jl")

Literate.markdown(
    INPUT,
    DOCS;
    execute = true,
    documenter = false,
    credit = false,
)

markdown = joinpath(DOCS, "integration_4body_b2ddKpi.md")
pdf = joinpath(DOCS, "integration_4body_b2ddKpi.pdf")

run(`pandoc $markdown --pdf-engine=lualatex -V monofont=Menlo -o $pdf`)
println("Wrote ", pdf)
