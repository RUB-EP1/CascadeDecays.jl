### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 6b2ff5d2-5b60-11f1-1680-094427183c9c
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	
	using CascadeDecays
	using HadronicLineshapes
	using ThreeBodyDecays: SpinParity, RecouplingLS, VertexFunction, @jp_str
end

# ╔═╡ 2c2294d8-b18f-4a79-985d-61d9230e08ff
system = let
	m_p, m_K, m_pi, m_Xic = 0.938272, 0.493677, 0.139570, 2.46771

	CascadeSystem(
    	SystemSpins(1, 0, 0; two_h0 = 1),
    	SystemMasses(m_p, m_K, m_pi; m0 = m_Xic))
end

# ╔═╡ 9f4c4169-8652-4e0f-b917-b05f1b88d3f0
# weak Ξc: known final parities, mother parity unknown
weak = add_parities(system, '+', '-', '-'; P0 = UndefinedParity)

# ╔═╡ a2eb7e7c-9ce0-4596-beaf-dc3c2595808e
topology_Λ = DecayTopology(((1, 2), 3))   # bracket: ((1,2),3)

# ╔═╡ 319cde3f-e813-4400-9fa2-7f8364b2c860
topology_Kx = DecayTopology((1, (2, 3)))   # bracket: (1,(2,3))

# ╔═╡ e6be26c0-dd29-4f09-ba1d-0c4b1f82b6e9
md"""
## Decay chain: Λ(1520) in pK   —  Ξc → [p K]_Λ π
"""

# ╔═╡ b805e622-03b6-46f9-8250-ace8caeeb65b
propagators_Λ = (
    (1, 2) => (
        jp = jp"1/2-",           # Λ(1520)−, J^P = 1/2−
        lineshape = BreitWigner(1.5195, 0.015),
    ),
)

# ╔═╡ 8eb77f47-2e51-4680-8288-bd924957a4fe
minimal_vertex_couplings(topology_Λ, weak, propagators_Λ)

# ╔═╡ 9b574a94-7de6-41ef-a8ed-f42c81a43d2f
chain_Λ  = minimal_ls_decay_chain(topology_Λ, weak, propagators_Λ)

# ╔═╡ 95a0ed60-396b-4f69-a3fe-347fd9ef59a8
[l.h.two_ls for l in chain_Λ.vertices]

# ╔═╡ bb81807e-65fd-4a78-ad61-abfa12cc2757
chains_Λ  = all_ls_decay_chains(topology_Λ, weak, propagators_Λ) |> size

# ╔═╡ f4d8f532-a916-4c88-9d5d-c26a2d35f090
md"""
## Option B: K*(892) in Kπ   —  `Ξc → p [K π]_K*`
"""

# ╔═╡ ea65535f-1a4d-4277-8825-916760f70716
propagators_Kx = (
    (2, 3) => (
        jp = jp"1-",           # K*(892)0, J^P = 1−
        lineshape = BreitWigner(0.89555, 0.047),
    ),
)

# ╔═╡ 46ec5dfd-fb4b-4d73-856a-3aed92ec5fb8
DecayChain(topology_Kx;
	propagators = propagators_Kx,
	vertices = (
		(1,(2,3)) => VertexFunction(RecouplingLS((2, 2))),
		(2,3) => VertexFunction(RecouplingLS((2, 2))),
	)
)

# ╔═╡ 237728d1-e197-4529-8538-caaa94c5325b
let
	(l1, itr_v1), (l2, itr_v2) = possible_vertex_couplings(topology_Kx, weak, propagators_Kx)
	map(Iterators.product(itr_v1, itr_v2)) do ((two_l1, two_s1), (two_l2, two_s2))
		# 
		DecayChain(topology_Kx;
			propagators = propagators_Kx,
			vertices = (
				l1 => VertexFunction(RecouplingLS((two_l1, two_s1))),
				l2 => VertexFunction(RecouplingLS((two_l1, two_s2))),
			)
		)
	end
end

# ╔═╡ f3e8fede-e12b-4d5d-9093-e338a8f541a4
chain_Kx = minimal_ls_decay_chain(topology_Kx, weak, propagators_Kx)

# ╔═╡ 7bcb0241-6226-4b0f-88cd-616f511754f2
chains_Kx = all_ls_decay_chains(topology_Kx, weak, propagators_Kx)

# ╔═╡ Cell order:
# ╠═6b2ff5d2-5b60-11f1-1680-094427183c9c
# ╠═2c2294d8-b18f-4a79-985d-61d9230e08ff
# ╠═9f4c4169-8652-4e0f-b917-b05f1b88d3f0
# ╠═a2eb7e7c-9ce0-4596-beaf-dc3c2595808e
# ╠═319cde3f-e813-4400-9fa2-7f8364b2c860
# ╟─e6be26c0-dd29-4f09-ba1d-0c4b1f82b6e9
# ╠═b805e622-03b6-46f9-8250-ace8caeeb65b
# ╠═8eb77f47-2e51-4680-8288-bd924957a4fe
# ╠═9b574a94-7de6-41ef-a8ed-f42c81a43d2f
# ╠═95a0ed60-396b-4f69-a3fe-347fd9ef59a8
# ╠═bb81807e-65fd-4a78-ad61-abfa12cc2757
# ╟─f4d8f532-a916-4c88-9d5d-c26a2d35f090
# ╠═ea65535f-1a4d-4277-8825-916760f70716
# ╠═46ec5dfd-fb4b-4d73-856a-3aed92ec5fb8
# ╠═237728d1-e197-4529-8538-caaa94c5325b
# ╠═f3e8fede-e12b-4d5d-9093-e338a8f541a4
# ╠═7bcb0241-6226-4b0f-88cd-616f511754f2
