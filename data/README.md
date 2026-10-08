# Exact-rational witness

`nongeneralizing_witness.txt` contains 135 rational coordinates for p=5, K=9, ordered neuron by neuron: alpha[1:5], beta[1:5], v[1:5]. The rationals are the exact dyadic values of a numerically discovered Float64 candidate.

The coordinates are an **approximate root**, not claimed to solve the equations exactly. Running

```sh
julia --project=. scripts/run.jl certify
```

proves the existence of an exact real interpolant in a small box around them. All certification inequalities use exact rational arithmetic. No binary Julia serialization is needed.

The observed pairs omit (2,4), (4,4), (5,4) in one-based array indices. The witness was regenerated during the October 9, 2026 reconstruction; the original temporary binary was lost. See `docs/PROVENANCE.md` and `docs/REPRODUCED.md`.
