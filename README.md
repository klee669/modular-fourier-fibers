# Fourier realizations and local completion in modular addition

This repository presents three computational results for the quadratic network

```text
Phi(alpha,beta,v)[a,b,c] = sum_j v[c,j] * (alpha[a,j] + beta[b,j])^2.
```

The target is the modular addition table. Array indices are one-based: the target is 1 when `mod(a+b-2,p)==c-1`, and 0 otherwise. Every observed input pair includes all p output channels; residuals are uncentered.

We study two explicit realizations of the same function: **eight neurons per conjugate frequency pair** (K=4p-3) and **four neurons per conjugate frequency pair** (K=2p-1). Each also has one constant neuron. These are descriptions of our constructions, not names of regimes from another paper. Their parameter-space geometries need not be the same.

## Run the three results

Install Julia **1.12.6**, then run from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run.jl
```

The default is `showcase`. 

```sh
julia --project=. scripts/run.jl showcase
```

The run prints a brief summary and writes exactly three report files:

```text
results/eight_neurons_per_frequency.toml
results/four_neurons_per_frequency.toml
results/nearby_completion.toml
```

Each report groups all of that experiment's cases into one TOML file. Running again overwrites these three files. It does not delete files produced by explicitly requested additional experiments. Set `MFF_RESULTS=/path/to/output` to use another output directory.

Nemo/FLINT install automatically. The exact characteristic-zero calculations can take some time; first use also includes Julia compilation.

## 1. Eight neurons per frequency pair: output-frequency blocks

Command for this result alone:

```sh
julia --project=. scripts/run.jl eight
```

Using `(u+v)^2-(u-v)^2=4uv`, the construction implements each Fourier product separately. It uses K=4p-3 neurons in total.

At this realization, exact characteristic-zero calculations over cyclotomic number fields show that the tangent image splits into independent output-frequency blocks:

| p | K | Constant output block dimension | Each nonconstant output block dimension | Full Jacobian rank |
|---:|---:|---:|---:|---:|
| 5 | 17 | 17 | 25 | 117 |
| 7 | 25 | 25 | 41 | 271 |

There are p-1 nonconstant real output modes. In each listed case, the sum of projected block dimensions equals the full rank. The patterns `4p-3` and `8p-15` are not yet proved for general p.

The report also restricts the p=7 blocks to partial observations. For the supplied seed-1 sample sets, 30 observed pairs give training rank 205, while the naive expression `min(7*30,271)` gives 210. The observed rank is the sum of the **restricted** block ranks. Thus both the number and placement of observations matter.

**Meaning:** recovering all first-order output variations requires sufficiently informative observations in every frequency block. This calculation does not establish nonlinear local completion for this construction. The listed ranks and direct-sum identities are exact in characteristic zero.

## 2. Four neurons per frequency pair: kernels and partial observations

Command for this result alone:

```sh
julia --project=. scripts/run.jl four
```

This construction combines Fourier products using identities such as

```text
(cos A + cos B)^2 + (sin A - sin B)^2 = 2 + 2 cos(A+B).
```

It uses K=2p-1 neurons. Its parameter directions couple output frequencies, so the independent-block result from the eight-neuron construction must not be transferred to it.

Exact characteristic-zero number-field computations give:

| p | K | Full rank | Kernel dimension | Translation/scaling gauge | Additional directions |
|---:|---:|---:|---:|---:|---:|
| 5 | 9 | 105 | 30 | 18 | 12 |
| 7 | 13 | 229 | 44 | 26 | 18 |

Each frequency pair supports four additional kernel directions beyond its eight translation/scaling directions. A further 4 directions at p=5, and 6 at p=7, lie outside the sum of the individual neuron-block kernels. These are infinitesimal statements, not a construction of continuous solution families.

For p=5, the report also compares the full Jacobian with the training Jacobian on two fixed observation sets:

- 22 pairs: omit `(2,4), (4,4), (5,4)` from the full table.
- 23 pairs: omit `(2,4), (5,4)`.

Both training ranks are 105, equal to the full rank. Since the full system contains every training equation, their parameter-space kernels coincide.

**Meaning:** at this Fourier point, every parameter direction preserving the observed outputs to first order also preserves the full table to first order. This does not itself prove smoothness or equality of the nonlinear fiber germs.

## 3. Nearby completion: numerical partial-fiber samples

Command for this result alone:

```sh
julia --project=. scripts/run.jl nearby
```

This experiment uses the **four-neuron-per-frequency-pair** construction at p=5, K=9, and the same two observation sets as Result 2. It:

1. Samples a unit direction in the training Jacobian kernel at the Fourier point.
2. Moves by a distance of 1e-2, 1e-3, or 1e-4.
3. Applies Newton correction using **only the observed-pair residuals**.
4. Measures both training and full-table residuals at the resulting point.

Two observation sets, three scales, and eight directions give 48 samples. The previous reproduction had maximum training and full residuals of about 1.95e-13. The new report records every residual and displacement rather than assuming this value will be identical on all platforms.

**Meaning:** these sampled nearby training solutions also fit the held-out entries to numerical precision. Together with Result 2 this motivates a local completion conjecture. Sampling is not a proof that every nearby training solution completes the table.

## Reading the results

Results 1 and 2 concern **different explicit realizations**. Result 3 is performed only around the realization in Result 2. We do not infer its nonlinear behavior from the independent-block decomposition in Result 1.

All three experiments concern completion at fixed p, not transfer to a new modulus or a larger group. The code was reconstructed from preserved exploratory records and rerun; the dependency graph is pinned in `Manifest.toml`.

The `src/` modules contain the computations and `scripts/run.jl` selects the experiments. Other research commands remain available only when explicitly requested; none are part of the default presentation. Older records in `docs/reference/` are archival and are not the current output directory.
