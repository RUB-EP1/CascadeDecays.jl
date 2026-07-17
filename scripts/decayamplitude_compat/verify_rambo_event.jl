using Random
using RamboOnDiet
using TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FIXTURE = joinpath(ROOT, "docs", "fixtures", "decayamplitude_compat", "event.toml")

event = TOML.parsefile(FIXTURE)
rng = MersenneTwister(event["seed"])
generator = PhaseSpaceGenerator(event["daughter_masses"], event["parent_mass"])
point = rand(rng, generator)
expected = sort(event["momenta"]; by = row -> row["particle"])

errors = Float64[]
for (momentum, row) in zip(point.momenta, expected)
    append!(errors, abs.(collect(momentum) .- [row["px"], row["py"], row["pz"], row["E"]]))
end
max_error = maximum(errors)

@assert max_error < 1.0e-14
println("RamboOnDiet fixture verification")
println("  seed: ", event["seed"])
println("  max four-vector error: ", max_error)
