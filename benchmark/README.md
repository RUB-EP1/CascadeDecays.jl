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

### Results (median time, 1000 events; both full matrices `(2,1,1,2)`)

| Chain | ThreeBodyDecays full | CascadeDecays full | Cascade / TBD (full) |
|-------|----------------------|--------------------|----------------------|
| K* | ~4.3 ms (~4.3 µs/event) | ~3.8 ms (~3.8 µs/event) | ~0.90× |
| Λ(1520) | ~5.0 ms (~5.0 µs/event) | ~4.0 ms (~4.0 µs/event) | ~0.82× |

Single helicity `(λ_p, λ_K, λ_π, λ_Λc) = (1,0,0,1)`: ThreeBodyDecays ~1.5–2.1 µs/event; Cascade ~3.8–4.1 µs/event (indexes the full array each call).

Both packages return **`(2, 1, 1, 2)`** on `(λ_p, λ_K, λ_π, λ_Λc)` once `ThreeBodySpins` uses `two_h0=1` (not `h0=1`, which incorrectly doubles every entry to spin 1). Internal resonance helicities are summed inside `ThreeBodyDecays.amplitude(dc, σs)` as in Cascade.

Single-helicity component uses `(λ_p, λ_K, λ_π, λ_Λc) = (1, 0, 0, 1)`.

After typing the internal contraction as `Val(Nf+1)`, `Val(Np)` on `DecayChain{Nf,Np}`, timings were unchanged within run-to-run noise (constant-folding already applied for concrete chains).
