# Lean verification status

Updated 2026-11-23. A statement appears in the first table only if it compiles in CI on the pinned toolchain (Lean 4.35.0, Mathlib commit 56a0058) with no `sorry` and with only the axioms `propext`, `Classical.choice` and `Quot.sound`.

## Verified in Lean

| Lean name | File | Statement | Supports |
| --- | --- | --- | --- |
| `scaledDet_mul_sin`, `scaledDet_eq`, `scaledDet_eq_zero_iff` | `Discrete/Recurrence.lean` | The recurrence $D_{n+2}=2\cos\theta\,D_{n+1}-D_n$, $D_0=1$, $D_1=2\cos\theta$ has $D_n\sin\theta=\sin((n+1)\theta)$; $D_n=0$ exactly when $\sin((n+1)\theta)=0$ | Determinant formula and discrete caustics, correspondence 002, Section 3 |
| `cos_sub_cos_ge` | `Discrete/SpectralGap.lean` | For $0\le a\le b$, $a+b\le\pi$: $\cos a-\cos b\ge2(b^2-a^2)/\pi^2$ | Spectral lower bound in Lemma D (002, Section 5) |
| `cos_gap_ge` | `Discrete/SpectralGap.lean` | For $2\le k$, $k+1\le N$: $\cos(\pi/N)-\cos(k\pi/N)\ge2(k^2-1)/N^2$ | Gap $\lambda_k-\lambda_1$ in the stiff Green bound (011, Section 3) |
| `sum_inv_sq_sub_one_eq`, `sum_inv_sq_sub_one_le` | `Discrete/SpectralGap.lean` | $\sum_{k=2}^{n+1}1/(k^2-1)=3/4-(2n+3)/(2(n+1)(n+2))\le3/4$ | Constant in $M\le3T/(8m)$ (011, Section 3) |
| `certificate_sum`, `certificate_two` | `Certificate.lean` | Logical assembly of the selection certificate from local bounds and a tail bound | Theorem S (008) |
| `first_term_alone_error_ge` | `Certificate.lean` | If the two-term formula has error at most $B$, the first term alone has error at least $\|L_2\|-B$ | Necessity statement (010) |

## Verified by interval arithmetic, not in Lean

| Item | Where | Remark |
| --- | --- | --- |
| Hypotheses of Theorem S at one parameter point, $\hbar=1/40,1/80,1/160$ | `experiments/2026-11-23-sparks/certify_polygon.py`, correspondence 010 | Arb ball arithmetic; program run by Sparks only |

## Conventional proof only

| Statement | Correspondence | Reviewed by |
| --- | --- | --- |
| Link between the recurrence and `Matrix.det` of the tridiagonal Hessian; eigenvalues and eigenvectors | 002 | Nova (003), numerics |
| Regulated Fresnel lemma | 002 | Nova (003) |
| Theorem C, Fourier expansion bound | 002 | Nova (003), numerics |
| Lemma D, uniqueness of the real critical point | 002 | Nova (003) |
| Connected-tree cumulant bound and zero-free region | 003 | Sparks (004), numerics |
| Proposition G, order $\hbar$ remainder | 004 | Nova (005) |
| Identification of the leading terms | 004, 005 | Open: one classical ingredient not written out |
| Proposition H, three real critical points | Nova, pull request 4 | Sparks, numerics |
| Lemma 1 (contour freedom), Lemma L (local contribution), Theorem S analytic inputs | 008 | Nova (009) |
| Three-saddle contour existence | 009 | Sparks (010) |
| Uniform soft-mode reduction | 011 | Sparks (013), numerics |
| Spectral truncation with pair correction | 012 | Sparks (013), not tested numerically |

## Next Lean targets, in order

1. The real-variable core of Lemma L.
2. The one-dimensional regulated Fresnel integral and its modulus bound.
3. Eigenvalues of the tridiagonal Hessian and the link to `Matrix.det`.
4. The moment sums of correspondence 011 and the Frobenius bound of correspondence 012, which are finite algebra.
