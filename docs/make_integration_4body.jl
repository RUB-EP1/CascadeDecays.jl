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

run(`pandoc $markdown --pdf-engine=lualatex -V monofont=Menlo -V geometry:left=2cm -V geometry:right=1.5cm -V geometry:top=1.5cm -V geometry:bottom=1.5cm -o $pdf`)
println("Wrote ", pdf)
