# Amplitude benchmarks

## Lc2pKpi (Λc⁺ → p K⁻ π⁺)

1000 weighted phase-space points from [RemboOnDiet.jl](https://github.com/mmikhasenko/RamboOnDiet.jl), masses from the package tutorial.

Decay chains use unit lineshape coupling (`σ -> 1`) and LS vertices matched to:

| Chain | ThreeBodyDecays | CascadeDecays topology |
|-------|-----------------|------------------------|
| K* in (p,K), π spectator | `DecayChainLS(k=3, jp="1-")` | `((1,2),3)` |
| Λ(1520) in (p,π), K spectator | `DecayChainLS(k=2, jp="3/2-")` | `((1,3),2)` |

Run:

```bash
julia --project=benchmark benchmark/Lc2pKpi.jl
```

### Results (median time, 1000 events)

| Chain | ThreeBody full | Cascade full | Ratio |
|-------|----------------|--------------|-------|
| K* | ~1.30 ms | ~2.65 ms | ~2.0× |
| Λ(1520) | ~1.85 ms | ~2.99 ms | ~1.6× |

Scalar helicity (all spin 0) timings are essentially the same as full matrix for this process.

After typing the internal contraction as `Val(Nf+1)`, `Val(Np)` on `DecayChain{Nf,Np}`, timings were unchanged within run-to-run noise (constant-folding already applied for concrete chains).
