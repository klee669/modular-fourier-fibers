# Conventions and certification

## Coordinates and realizations

Array indices a,b,c run from 1 through p, representing group elements a-1,b-1,c-1. Thus the target entry equals 1 exactly when `mod(a+b-2,p)==c-1`.

Parameters are neuron-major: for each neuron j, concatenate `alpha[:,j]`, `beta[:,j]`, `v[:,j]`. Outputs are channel-major, then ordered pairs `(1,1),(1,2),...,(p,p)`, or the specified ordering of D. No gauge fixing or centering is implicit in the Jacobian.

For each positive frequency k, set C(a)=cos(2 pi k a/p), S(a)=sin(2 pi k a/p). The narrow construction uses the identities

```text
(C(a)+C(b))^2 + (S(a)-S(b))^2 = 2 + 2 cos(2 pi k(a+b)/p)
(S(a)+C(b))^2 + (C(a)+S(b))^2 = 2 + 2 sin(2 pi k(a+b)/p).
```

Their output coefficients are C(c)/p and S(c)/p. One constant neuron corrects the constant main effects and supplies the constant Fourier mode. There are `4*(p-1)/2+1=2p-1` neurons. The wide realization instead uses `(a+b)^2-(a-b)^2=4ab` for each of the four Fourier products, with `4p-3` neurons.

The reconstruction consistently places the constant neuron last. The original wide scripts sometimes placed it first. This changes parameter column order, not ranks or tangent images.

Translation gauge: `delta alpha[:,j]=1`, `delta beta[:,j]=-1`, `delta v[:,j]=0`.
Scaling gauge: `delta alpha[:,j]=alpha[:,j]`, `delta beta[:,j]=beta[:,j]`, `delta v[:,j]=-2v[:,j]`.

## Obstruction computation

At a root x0 write the cubic residual as

```text
F(x0+h) = J h + Q(h,h) + C(h,h,h).
Q(h,h) = sum_j [2 dv_j z_j dz_j + v_j dz_j^2].
C(h,h,h) = sum_j dv_j dz_j^2.
```

Here z is the **base** preactivation alpha_a+beta_b. Its factor in the mixed quadratic term is essential. Tests check the exact Taylor identity at a random base point and check coefficient tables by substitution.

Let E span a complement of translation/scaling gauges in the common Jacobian kernel. Let L have columns spanning the left kernel of J. The quadratic obstruction is `L' Q(E*s,E*s)`. If it vanishes, solve coefficientwise for a quadratic correction Y(s) with

```text
J Y(s) = -Q(E*s,E*s).
```

The corrected cubic obstruction is then

```text
L' [2 Q(E*s,Y(s)) + C(E*s,E*s,E*s)].
```

The code includes all off-diagonal monomial multiplicities. For a 12-dimensional slice there are 78 degree-two and 364 degree-three monomials. Partial-data corrections are also checked against full-data equations at degree two. The tests do not claim vanishing in all orders.

## Exact existence and submersion certificate

Take a numerical candidate x, convert its Float64 coordinates exactly to dyadic rationals, and select m=p|D| parameter coordinates. Fix every other coordinate at its rational value. This defines a square map f:R^m -> R^m. Let B be the **exact rational inverse** of the selected Jacobian at x.

Use the maximum norm and its induced matrix norm. Set

```text
eta = ||B f(x)||_infinity
rho = 2 eta
M   = max_i |x_i| + rho
L   = 24 K M
q   = ||B||_infinity L rho.
```

For each scalar output of one neuron, the sum of absolute Hessian entries is at most

```text
8 |v| + 8 |z| <= 24 M,  since |z| <= 2M.
```

Summing over neurons and restricting to selected coordinates gives the Jacobian Lipschitz bound L throughout the closed rho-box. All entries of x, B, eta, rho, M, L, q are rational; the comparisons are exact.

If `q<1` and `eta+q*rho<rho`, then the map `G(y)=y-B f(y)` maps the box into itself and is a contraction. Banach's theorem gives a unique zero within the **selected affine slice** and box. The full parameter-space fiber is not claimed to be a single point.

Moreover, `||I-B Df(y)|| <= q < 1` throughout the box. Hence the selected square derivative is invertible everywhere in it. The original training map therefore has full row rank at the certified zero, and its fiber is locally smooth of dimension `3pK-p|D|`.

For a held-out output the sum of absolute first derivatives is bounded by

```text
sum_j (4 |v_j z_j| + |z_j|^2) <= 12 K M^2.
```

Consequently, `|F_heldout(x)| - 12 K M^2 rho > 0` certifies a nonzero held-out error at the exact root. No floating-point residual tolerance is used as an existence proof.

For uncentered sum-of-squares loss, submersion and a smooth positive local prior give local learning coefficient m/2. No general theorem identifying LLC with half the Jacobian rank at a singular point is assumed.
