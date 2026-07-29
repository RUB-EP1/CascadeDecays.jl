### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ c47c14b0-0c65-4a50-b69d-f92b48fc5e3f
begin
    using Pkg
    Pkg.activate(@__DIR__)

    using HypertextLiteral: @htl, JavaScript
    using JSON
    using PlutoUI

    include(joinpath(@__DIR__, "lb2lc3pi_vortex_strings.jl"))
end

# ╔═╡ d0e1292d-b68d-4f33-8ef8-7144c461f46f
md"""
# Vortex strings in $\Lambda_b^0\to\Lambda_c^+\pi^+\pi^-\pi^-$

For the two-by-two external-helicity amplitude matrix $A$, the condition

```math
\det A = 0
```

is two real equations: $\operatorname{Re}\det A=0$ and
$\operatorname{Im}\det A=0$. The orientation-free four-body phase space has
five dimensions, so a regular zero has dimension $5-2=3$. Fixing the two
relative angles leaves a three-mass volume; the same zero set then appears as
**one-dimensional vortex strings**.
"""

# ╔═╡ da647431-d3f1-4475-8f0a-c02850af57d8
md"""
## A topology-adapted square chart

Choose a sequential coordinate tree

```text
Λb → (a,(b,c)) + d .
```

The displayed coordinates are the three masses
$(m_{ab},m_{ac},m_{abc})$. Internally, every point is generated from five
bounded variables:

```text
u_cluster, t_pair, v_helicity, u_theta, u_phi ∈ (0,1).
```

They map respectively to $m_{abc}$, $m_{bc}$, $\cos\chi$, $\theta_a$, and
$\phi_b$. Thus every topology has the same square computational domain even
though its physically natural masses are different.
"""

# ╔═╡ b136814f-f77f-4fd6-88f0-51209445f3ba
chart_options = [String(chart.key) => chart.bracket for chart in VORTEX_CHARTS]

# ╔═╡ 3740fb74-f4a8-4ceb-bff3-fd50cb988464
@bind selected_chart_key Select(chart_options; default = "a1_34")

# ╔═╡ e8b93971-835f-4fc3-bddb-ab9d846fbba8
selected_chart = vortex_chart(selected_chart_key)

# ╔═╡ 8c30d119-63e7-4bbc-b94a-194311991597
chart_info = chart_description(selected_chart)

# ╔═╡ 191c69e7-40dd-492b-a91d-315f21bf8d06
md"""
Selected topology: **$(chart_info.topology)**

- masses: `$(join(chart_info.mass_labels, ", "))`
- angles: `$(join(chart_info.angle_labels, ", "))`
- bachelor: $(chart_info.particle_labels.bachelor)
- pair: $(join(chart_info.particle_labels.pair, ", "))
- spectator: $(chart_info.particle_labels.spectator)
"""

# ╔═╡ 1e998041-58fa-47a2-8c6e-ae0904bdab75
@bind controls confirm(
    PlutoUI.combine() do Child
        @htl("""
        <div style="display:grid; grid-template-columns:11rem minmax(18rem,1fr);
                    gap:.55rem 1rem; align-items:center; max-width:46rem">
          <label>cluster mass, u</label>
          $(Child(:u_cluster, Slider(0.03:0.01:0.97; default=0.93, show_value=true)))
          <label>pair mass, t</label>
          $(Child(:t_pair, Slider(0.02:0.01:0.98; default=0.50, show_value=true)))
          <label>helicity, v</label>
          $(Child(:v_helicity, Slider(0.02:0.01:0.98; default=0.50, show_value=true)))
          <label>polar angle, uθ</label>
          $(Child(:u_theta, Slider(0.02:0.01:0.98; default=1/3, show_value=true)))
          <label>azimuth, uϕ</label>
          $(Child(:u_phi, Slider(0.02:0.01:0.98; default=0.75, show_value=true)))
          <label>mass-window fraction</label>
          $(Child(:mass_span, Slider(0.08:0.02:0.50; default=0.14, show_value=true)))
          <label>Dalitz planes</label>
          $(Child(:nslices, Slider(3:2:15; default=7, show_value=true)))
          <label>initial root grid</label>
          $(Child(:ngrid, Slider(5:2:11; default=7, show_value=true)))
        </div>
        """)
    end;
    label = "Trace vortex strings",
)

# ╔═╡ 5ee21fb6-6b2a-45d4-8fc2-72c8dfde77f7
coordinates = let
    limits = chart_limits(selected_chart)
    mabc = limits.lower + controls.u_cluster * (limits.upper - limits.lower)
    pair_fraction = controls.t_pair
    cos_helicity = 2controls.v_helicity - 1
    theta = π * controls.u_theta
    phi = π * (2controls.u_phi - 1)
    masses = chart_masses(
        selected_chart, mabc, pair_fraction, cos_helicity,
    )
    (; limits, mabc, pair_fraction, cos_helicity, theta, phi, masses)
end

# ╔═╡ 0e2e5d56-375b-48c2-939e-840299f73852
model, task = lb2lc3pi_model()

# ╔═╡ 7d11b47c-72a8-4737-823e-48c2e4be20fd
probe_det = chart_normalized_det(
    model,
    task,
    selected_chart,
    coordinates.mabc,
    coordinates.pair_fraction,
    coordinates.cos_helicity,
    coordinates.theta,
    coordinates.phi,
)

# ╔═╡ 5e16f653-e9f0-4d8a-bdde-dde8af5eed4d
md"""
### Current phase-space point

| coordinate | value |
|:--|--:|
| $(chart_info.mass_labels[1]) | $(round(coordinates.masses.mab; digits=6)) GeV |
| $(chart_info.mass_labels[2]) | $(round(coordinates.masses.mac; digits=6)) GeV |
| $(chart_info.mass_labels[3]) | $(round(coordinates.mabc; digits=6)) GeV |
| $(chart_info.angle_labels[1]) | $(round(rad2deg(coordinates.theta); digits=2))° |
| $(chart_info.angle_labels[2]) | $(round(rad2deg(coordinates.phi); digits=2))° |
| ``\left|\det A/\sum |A|^2\right|`` | $(round(abs(probe_det); sigdigits=6)) |
"""

# ╔═╡ 4d5d5bff-771d-4133-b774-a433d518441d
window = let
    full_width = coordinates.limits.upper - coordinates.limits.lower
    half_width = controls.mass_span * full_width / 2
    epsilon = 2e-3
    mlo = max(coordinates.limits.lower + epsilon, coordinates.mabc - half_width)
    mhi = min(coordinates.limits.upper - epsilon, coordinates.mabc + half_width)
    (; mlo, mhi)
end

# ╔═╡ f1020b04-fb6f-48b0-8864-129d13524f65
vortices = trace_chart_strings(
    model,
    task,
    selected_chart;
    theta = coordinates.theta,
    phi = coordinates.phi,
    nslices = Int(controls.nslices),
    ngrid = Int(controls.ngrid),
    continuation_ngrid = 3,
    window.mlo,
    window.mhi,
)

# ╔═╡ 23832bc9-85f5-4fb2-a64c-421408399d0c
plot_data = let
    borders = [
        chart_dalitz_border(selected_chart, mass; n = 90)
        for mass in vortices.mgrid
    ]
    branch_numbers = sort(unique(p.branch for p in vortices.points))
    Dict(
        "labels" => collect(chart_info.mass_labels),
        "borders" => [
            Dict(
                "x" => [p.mab for p in border],
                "y" => [p.mac for p in border],
                "z" => fill(vortices.mgrid[i], length(border)),
            )
            for (i, border) in enumerate(borders)
        ],
        "vortices" => [
            let points = filter(p -> p.branch == branch, vortices.points)
                Dict(
                    "branch" => branch,
                    "x" => [p.mab for p in points],
                    "y" => [p.mac for p in points],
                    "z" => [p.mabc for p in points],
                    "residual" => [p.residual for p in points],
                )
            end
            for branch in branch_numbers
        ],
        "probe" => Dict(
            "x" => [coordinates.masses.mab],
            "y" => [coordinates.masses.mac],
            "z" => [coordinates.mabc],
        ),
    )
end

# ╔═╡ a5bbfa81-018a-402b-85aa-59fdd1ef729a
let
    payload = JavaScript(JSON.json(plot_data))
    @htl("""
    <div style="border:1px solid #d7dce2; border-radius:12px; padding:.4rem;
                background:#fbfcfe">
      <div class="vortex-plot" style="width:100%; height:680px"></div>
      <script>
        const payload = $payload
        const plotNode = currentScript.parentElement.querySelector(".vortex-plot")
        const module = await import(
          "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@3.0.1/+esm"
        )
        const Plotly = module.default ?? module
        const borderTraces = payload.borders.map((border, index) => ({
          type: "scatter3d",
          mode: "lines",
          x: border.x, y: border.y, z: border.z,
          line: {color: "rgba(55,74,96,.56)", width: 3},
          name: index === 0 ? "Dalitz borders" : undefined,
          showlegend: index === 0,
          hoverinfo: "skip"
        }))
        const vortexTraces = payload.vortices.map((branch, index) => ({
          type: "scatter3d",
          mode: "markers+lines",
          x: branch.x,
          y: branch.y,
          z: branch.z,
          marker: {color: "#e32636", size: 5},
          line: {color: "#e32636", width: 3},
          customdata: branch.residual,
          hovertemplate:
            payload.labels[0] + "=%{x:.5f}<br>" +
            payload.labels[1] + "=%{y:.5f}<br>" +
            payload.labels[2] + "=%{z:.5f}<br>" +
            "|normalized det|=%{customdata:.2e}<extra>vortex</extra>",
          name: "det(A) = 0",
          showlegend: index === 0
        }))
        const probeTrace = {
          type: "scatter3d",
          mode: "markers",
          x: payload.probe.x, y: payload.probe.y, z: payload.probe.z,
          marker: {color: "#f2a900", size: 8, symbol: "diamond"},
          name: "slider point"
        }
        Plotly.react(plotNode, [...borderTraces, ...vortexTraces, probeTrace], {
          margin: {l: 0, r: 0, b: 0, t: 25},
          paper_bgcolor: "#fbfcfe",
          scene: {
            xaxis: {title: payload.labels[0] + " [GeV]"},
            yaxis: {title: payload.labels[1] + " [GeV]"},
            zaxis: {title: payload.labels[2] + " [GeV]"},
            aspectmode: "data",
            camera: {eye: {x: 1.55, y: 1.45, z: 1.15}}
          },
          legend: {x: .01, y: .99}
        }, {responsive: true, displaylogo: false})
      </script>
    </div>
    """)
end

# ╔═╡ f0b4ec16-a68a-4499-a4db-1702a7db7d1e
md"""
The gray loops are physical Dalitz boundaries at several values of the cluster
mass. Red points are numerical solutions of
$\det A/\sum_{\lambda',\lambda}|A_{\lambda'\lambda}|^2=0$; red segments join
solutions found by continuation between neighboring mass slices. The gold
diamond is the independent point controlled by all five sliders.

Found **$(length(vortices.points))** vortex intersections on
**$(length(vortices.mgrid))** slices. The largest root residual is
**$(isempty(vortices.points) ? "n/a" :
    string(round(maximum(p.residual for p in vortices.points); sigdigits=3)))**.
"""

# ╔═╡ a4fbd91c-6241-4893-b370-78fdff02f28b
md"""
## How the calculation is performed

1. In the $\Lambda_b$ rest frame, the $(abc)$ cluster is placed along $+z$.
2. In the cluster rest frame, bachelor $a$ is put in the $xz$ plane at
   $\theta_a$ with positive azimuth convention.
3. In the $(bc)$ rest frame, particle $b$ is placed at $(\chi,\phi_b)$.
   Two boosts and two $y$ rotations reconstruct the four final momenta.
4. A fixed global rotation moves the event away from a helicity-coordinate
   pole. It changes neither invariants nor relative angles; without it, an
   exactly axial representative can create false zeros from undefined
   azimuths.
5. For each fixed $m_{abc}$ plane, a damped two-dimensional Newton solver
   finds simultaneous zeros of the real and imaginary determinant. Roots from
   one plane seed the next plane, which is the numerical parametrization of
   the string.

The normalization by $\sum|A|^2$ changes the scale but not any interior zero.
Near a regular root the real $2\times2$ Jacobian in the Dalitz plane has full
rank, and the determinant phase winds by an integer around a small loop—the
local vortex diagnostic.
"""

# ╔═╡ Cell order:
# ╠═c47c14b0-0c65-4a50-b69d-f92b48fc5e3f
# ╟─d0e1292d-b68d-4f33-8ef8-7144c461f46f
# ╟─da647431-d3f1-4475-8f0a-c02850af57d8
# ╠═b136814f-f77f-4fd6-88f0-51209445f3ba
# ╠═3740fb74-f4a8-4ceb-bff3-fd50cb988464
# ╠═e8b93971-835f-4fc3-bddb-ab9d846fbba8
# ╠═8c30d119-63e7-4bbc-b94a-194311991597
# ╟─191c69e7-40dd-492b-a91d-315f21bf8d06
# ╠═1e998041-58fa-47a2-8c6e-ae0904bdab75
# ╠═5ee21fb6-6b2a-45d4-8fc2-72c8dfde77f7
# ╠═0e2e5d56-375b-48c2-939e-840299f73852
# ╠═7d11b47c-72a8-4737-823e-48c2e4be20fd
# ╟─5e16f653-e9f0-4d8a-bdde-dde8af5eed4d
# ╠═4d5d5bff-771d-4133-b774-a433d518441d
# ╠═f1020b04-fb6f-48b0-8864-129d13524f65
# ╠═23832bc9-85f5-4fb2-a64c-421408399d0c
# ╠═a5bbfa81-018a-402b-85aa-59fdd1ef729a
# ╟─f0b4ec16-a68a-4499-a4db-1702a7db7d1e
# ╟─a4fbd91c-6241-4893-b370-78fdff02f28b
