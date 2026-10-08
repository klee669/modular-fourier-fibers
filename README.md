# Modular addition: Fourier fibers and certified interpolants

A standalone, reproducible computational study of

```text
Phi(alpha,beta,v)[a,b,c] = sum_j v[c,j] * (alpha[a,j] + beta[b,j])^2
T[a,b,c] = 1[(a+b) mod p = c].
```

This repository reconstructs exploratory calculations originally performed in September 2026. The original temporary directory was lost. The reconstruction includes new execution logs, meaningful regression tests, and a portable exact-rational witness. It is **not** a byte-for-byte recovery of the old scripts or a proof of general formulas.

## Reproduce

Install **Julia 1.12.6** (the tested version). Nemo/FLINT are installed automatically. No Python, Singular, GPU, or other local repository is needed.

From this repository's root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
julia --project=. scripts/run.jl quick
```

`quick` runs finite-field rank/block experiments, the moderate-data comparison, all quadratic/cubic obstruction coefficients on a gauge slice, phase/splitting checks, nearby numerical sampling, and exact certification of the supplied witness. Reports are written as machine-readable TOML under `results/`. First use includes compilation time.

```sh
# Add characteristic-zero number-field ranks at p=5,7 and larger block checks.
julia --project=. scripts/run.jl full

# Certify a supplied witness, without numerical optimization.
julia --project=. scripts/run.jl certify

# Rediscover a witness from random initializations, then certify it.
julia --project=. scripts/run.jl search
julia --project=. scripts/run.jl certify results/discovered_witness.txt
```

The search uses seed 1729 and one BLAS thread. Exact iterates can vary across Julia/BLAS versions. Certification of the bundled rational witness is independent of the search. Floating QR only selects a square minor: its inverse and all certification inequalities are subsequently computed exactly over Q.

Other commands: `ranks`, `moderate`, `jets`, `realizations`, `nearby`. Set `MFF_RESULTS=/path/to/output` to change the output directory. GitHub Actions runs tests and `quick` on Linux; this workflow is supplied, not a claim that remote CI has already run.

## What is being reproduced?

| Calculation | Report / command | Evidence |
|---|---|---|
| Moderate-data Fourier and nongeneralizing interpolants both have full-row-rank training Jacobians | `moderate` | Exact interpolation; rank bounds certified by finite-field minors |
| Wide Fourier tangent image splits by output frequency; naive rank saturation can fail | `ranks`, `wide_sweep_*.toml` | Exact finite-field linear algebra; finite tested cases |
| Narrow Fourier kernel = gauge plus additional directions | `kernel_*.toml` | Exact finite-field and, in `full`, characteristic-zero number-field ranks |
| Training/full tangent kernels agree for fixed p=5 observation sets | `jets` | Exact finite-field linear algebra |
| Quadratic and corrected cubic obstructions vanish on a 12-dimensional gauge slice | `jets` | All 78 quadratic and 364 cubic coefficients, not random-direction tests |
| Phase changes and neuron splitting | `realizations` | Floating phase tests; exact finite-field splitting test |
| Nearby partial-fiber Newton samples also fit held-out outputs | `nearby` | Numerical evidence only |
| A nongeneralizing submersive interpolant exists for the same p=5, K=9, D | `certify` | Exact-rational contraction certificate |

The fixed comparison problem has p=5, K=9, and 22 observed pairs. In **one-based** array coordinates it omits `(2,4), (4,4), (5,4)` from the full 5-by-5 table. The 23-pair variant omits only `(2,4), (5,4)`.

## Main reference values

| Realization | p | K | Full Jacobian rank | Kernel | Translation/scaling gauge |
|---|---:|---:|---:|---:|---:|
| narrow | 5 | 9 | 105 | 30 | 18 |
| narrow | 7 | 13 | 229 | 44 | 26 |
| wide | 5 | 17 | 117 | 138 | 34 |
| wide | 7 | 25 | 271 | 254 | 50 |

For the wide realization, the projected block dimensions at p=5,7,11 are respectively `[17,25,...,25]`, `[25,41,...,41]`, and `[41,73,...,73]`. The corresponding full ranks are 117,271,771. The proposed formulas `d0=4p-3`, `d1=8p-15` are **conjectural beyond verified cases** and require handling small-p collisions. Narrow projected blocks are coupled; adding their dimensions gives the wrong full rank.

For narrow p=7, the kernel decomposes into 26 gauge directions and 18 extra directions: 4 extra directions per each of three frequency blocks, plus 6 directions outside the sum of individual block kernels. This is an infinitesimal statement, not an integration theorem.

The certified nongeneralizing interpolant has training rank 110 and a smooth 25-dimensional training fiber. Under the uncentered squared-residual loss and a smooth positive local prior, its local learning coefficient is 55. A value 52.5 for the Fourier point remains **conditional** on proving the relevant 30-dimensional fiber is smooth with the expected tangent space.

## Scope and limitations

- Every observed pair includes all p output channels. Losses here are **uncentered**. These are not automatically results for sample-centered logits or regularized objectives.
- Fixed-p table completion is not transfer to a new modulus, a larger input size, or a new group.
- Finite-field computations are exact in their stated characteristic. They do not by themselves prove characteristic-zero upper bounds or identities.
- Vanishing obstructions through degree three and nearby Newton samples do not establish germ equality, smoothness, irreducibility, or RLCT at the Fourier point.
- The existence certificate concerns one nongeneralizing solution. It does not prove that it belongs to a different irreducible component, or that generalizing points globally minimize LLC.
- The cubic residual map can acquire higher-order terms after implicit elimination; testing degree three is not testing all orders.
- Permutations and finite cyclic group actions do not contribute infinitesimal gauge directions.

See [the mathematical conventions and certificate proof](docs/MATHEMATICS.md), [reconstruction provenance](docs/PROVENANCE.md), and [recorded reproduction results](docs/REPRODUCED.md).

## Layout

```text
src/                 model, exact linear algebra, jets, numerical search, certificate
scripts/run.jl       command-line experiments
test/runtests.jl     exact Taylor/polarization and monomial-coefficient checks; certificate
data/                portable witness (exact rational coordinates)
docs/                definitions, certificate proof, provenance, reproduction results
results/             generated TOML reports (ignored by git)
Manifest.toml        pinned dependency graph
```

The default runs do not modify the bundled witness. To contribute another witness, first run `search`, certify the generated text file, and record the associated environment and bounds.
