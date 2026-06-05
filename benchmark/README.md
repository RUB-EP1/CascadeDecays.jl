# Lc2pKpi dev script

Local smoke check for Λc⁺ → p K⁻ π⁺ against [ThreeBodyDecays.jl](https://github.com/gw2ssi/ThreeBodyDecays.jl) on the **same** K* and Λ(1520) LS chains. This is **not** a performance benchmark and **timing numbers must not be quoted** elsewhere.

## What is compared

- 1000 phase-space points from [RemboOnDiet.jl](https://github.com/mmikhasenko/RamboOnDiet.jl)
- Unit lineshape (`σ -> 1`), minimal LS vertices
- External helicity array shape **`(2, 1, 1, 2)`** for `(λ_p, λ_K, λ_π, λ_Λc)` when `ThreeBodySpins(1, 0, 0; two_h0 = 1)` (not `h0 = 1`)

| Chain | ThreeBodyDecays | CascadeDecays topology |
|-------|-----------------|------------------------|
| K* in (p,K), π spectator | `DecayChainLS(k=3, jp="1-")` | `((1,2),3)` |
| Λ(1520) in (p,π), K spectator | `DecayChainLS(k=2, jp="3/2-")` | `((1,3),2)` |

## Why timings are not comparable

`ThreeBodyDecays.amplitude(dc, σs)` returns the **helicity-frame** amplitude: it builds the aligned chain (`aligned_amplitude`, vertex × Wigner × vertex) and then applies **four Wigner `d` rotations** (parent + three finals) via `@tullio`.

`CascadeDecays.amplitude` stops after the **line buffer** and **internal helicity contraction** (aligned-style vertex product, no extra frame rotations). It does **not** include those helicity-frame Wigner steps yet. Any side-by-side `@benchmark` output compares different physics workloads; a smaller elapsed time for Cascade is **misleading**, not a speed win.

Use the script only to sanity-check shapes (e.g. `(2, 1, 1, 2)`, `A[:,1,1,:]` is `2×2` for the Λc–p block) while developing the cascade evaluator.

Run:

```bash
julia --project=benchmark benchmark/Lc2pKpi.jl
```
