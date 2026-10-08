# Reconstruction provenance

The September 2026 exploratory scripts lived in `/private/tmp/modular_fiber_probe.QY3buF`. On October 9 this directory was absent. This project was reconstructed from code excerpts and actual tool outputs preserved in the conversation, then rerun in a new standalone environment.

## Preserved versus regenerated

- Model, narrow/wide Fourier identities, fixed p=5 observation sets, gauge vectors, and obstruction formulas follow the surviving scripts.
- The p=5, K=9, n=22 certified comparison uses the same observation set: omit one-based pairs (2,4), (4,4), (5,4).
- The old serialized floating witness was lost. The bundled text witness is a **newly generated and independently certified** artifact, not a claim to recover the old binary. The search algorithm and seed follow the recorded successful search.
- Generic observation sweeps are regenerated using explicit seeds; they are not claimed to reproduce each old shuffled set. Every conclusion must be read from the newly generated reports.
- The original wide construction's constant neuron was sometimes first; the new code puts it last.
- The old finite-field "partial obstruction rank 4/full rank 8" result was caused by a missing base preactivation z in Q. It was corrected during the original work and is excluded from this project. Its downstream cone/radical computations are also excluded.
- An early conjecture `rank J_D = min(p|D|, rank J_all)` was contradicted by the wide block experiments. It is not asserted here.
- A conjectural generic rank formula for arbitrary weights is not established by these Fourier-point computations and is not presented as a result.

## What has not been proved

General-p formulas; smoothness or local dimension of the Fourier fiber; equality of partial and full fiber germs; integration of the additional kernel directions; irreducible-component membership; global minimum of LLC; and a causal account of grokking dynamics remain outside these computations.

Reports distinguish arithmetic over finite fields, arithmetic over exact number fields, numerical rank/sampling, and exact rational existence certificates. A finite-field equality is not silently promoted to a characteristic-zero theorem.

No remote repository has to be configured to run this project. It has no dependencies on the parent CertifiedHomotopyTracking checkout.
