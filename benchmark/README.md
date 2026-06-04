# Amplitude benchmarks

## Lc2pKpi (Λc⁺ → p K⁻ π⁺)

1000 weighted phase-space points from [RemboOnDiet.jl](https://github.com/mmikhasenko/RamboOnDiet.jl), masses from the package tutorial.

Hadronic quantum numbers: **p 1/2⁺**, **K 0⁻**, **π 0⁻**, **Λc 1/2⁺** (Λc parity fixed to `+` for LS chain building; use `jp0="1/2±"` in production if both parities are needed).

External helicity array shape for Cascade: **(2, 1, 1, 2)** — the **Λc–p block** `A[:,1,1,:]` is **2×2**.

Decay chains use unit lineshape coupling (`σ -> 1`) and minimal LS vertices (`minimal_ls_decay_chain`):

| Chain | ThreeBodyDecays | CascadeDecays topology |
|-------|-----------------|------------------------|
| K* in (p,K), π spectator | `DecayChainLS(k=3, jp="1-")` | `((1,2),3)` |
| Λ(1520) in (p,π), K spectator | `DecayChainLS(k=2, jp="3/2-")` | `((1,3),2)` |

Run:

```bash
julia --project=benchmark benchmark/Lc2pKpi.jl
```

### Results (median time, 1000 events)

| Chain | ThreeBody full `(3,1,1,3)` | Cascade full `(2,1,1,2)` | Ratio (Cascade/TBD) |
|-------|---------------------------|------------------------|---------------------|
| K* | ~17 ms | ~3.7 ms | ~0.2× |
| Λ(1520) | ~18 ms | ~4.3 ms | ~0.2× |

ThreeBodyDecays includes the resonance spin in the aligned tensor (larger `(3,1,1,3)`); Cascade keeps only external lines in the final array. Single-helicity component uses `(λ_p, λ_K, λ_π, λ_Λc) = (1, 0, 0, 1)`.

After typing the internal contraction as `Val(Nf+1)`, `Val(Np)` on `DecayChain{Nf,Np}`, timings were unchanged within run-to-run noise (constant-folding already applied for concrete chains).
