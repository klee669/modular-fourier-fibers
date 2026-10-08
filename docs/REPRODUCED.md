# Reproduction record — October 9, 2026

Environment: Julia 1.12.6, Nemo 0.56.2, AbstractAlgebra 0.50.3, macOS/aarch64. The exact dependency graph is in `Manifest.toml`. Numerical experiments use one BLAS thread.

Completed locally:

```sh
julia --project=. scripts/run.jl search
julia --project=. test/runtests.jl
julia --project=. scripts/run.jl full
```

All 28 test assertions passed: 12 root/gauge checks, 4 exact Taylor/polarization/coefficient-assembly checks, 7 rank/jet checks, and 5 portable-witness/certificate checks. The complete `full` command exited successfully. GitHub Actions is configured separately and has not been run remotely here.

## Exact ranks

Characteristic-zero number-field computations gave:

| p | K | full rank | kernel | gauge | extra | block-supported kernel dimensions | remaining coupled kernel |
|---:|---:|---:|---:|---:|---:|---|---:|
| 5 | 9 | 105 | 30 | 18 | 12 | 12,12,2 | 4 |
| 7 | 13 | 229 | 44 | 26 | 18 | 12,12,12,2 | 6 |

Both were also reproduced over two finite fields per p. For p=7 the primes were 1000133 and 2000153.

Finite-field full ranks for wide realizations at p=5,7,11 were 117,271,771. The reconstructed p=11 narrow realization had rank 621 in the tested finite field; this is not a characteristic-zero certification of that upper bound.

For a regenerated seed-1 nested-style size sweep (each observation set is specified by the public sampler), the wide p=7 ranks were:

| observed pairs | actual rank | min(7n,271) |
|---:|---:|---:|
| 20 | 140 | 140 |
| 25 | 173 | 175 |
| 30 | 205 | 210 |
| 35 | 235 | 245 |
| 40 | 265 | 271 |
| 45 | 271 | 271 |

At n=25 the constant output block itself has rank 23 rather than 25. Thus even the blockwise maximal-rank formula needs a sampling/unisolvence condition; it is not a formula for every D.

## Jets and numerical realizations

For each of the fixed p=5, K=9 observation sets of sizes 22 and 23, and each of two finite fields:

- training and full ranks were 105;
- the gauge complement had dimension 12;
- all 78 quadratic obstruction coefficients vanished;
- the quadratic correction from partial equations also solved full equations;
- all 364 corrected cubic obstruction coefficients vanished.

Three continuous phase variants had full and 22-pair ranks 105. Exact splitting of one neuron gave K=10, rank 105, kernel 45.

Across 48 nearby numerical projections (two D sets, three perturbation scales, eight draws), the largest observed training and full residual infinity norms were each approximately 1.954e-13. This is sampling evidence, not a local completion proof.

## Independent rediscovery and rational certification

The second initialization of the seed-1729 search reproduced the formerly recorded numerical candidate values:

- 23 Levenberg–Marquardt iterations;
- training residual infinity norm 9.345857421294568e-13;
- full residual infinity norm 10.193416419629049;
- smallest training Jacobian singular value 0.020926533993406292.

The new text witness is `data/nongeneralizing_witness.txt`. The exact rational certificate yielded:

| Bound or conclusion | Value |
|---|---:|
| existence-box radius | 2.0104907686938707e-10 |
| contraction bound | 2.6476347318071285e-5 |
| training rank at certified root | 110 |
| local training-fiber dimension | 25 |
| held-out coordinate, one-based | (5,4,4) |
| held-out absolute error lower bound | 10.19341636312218 |
| LLC under the stated loss/prior convention | 55 |

Decimal values above are displays. The certificate uses exact rational values and strict rational comparisons, included in [reference/certificate.toml](reference/certificate.toml). Its proof is in [MATHEMATICS.md](MATHEMATICS.md).

All generated TOML reports from this reconstruction are copied to `docs/reference/`. These are reference artifacts; fresh runs write to `results/` so they do not overwrite the record. Numerical values and selected QR columns may vary across platforms. Exact pass/fail certification and mathematical rank conclusions should not depend on matching floating-point output byte-for-byte.
