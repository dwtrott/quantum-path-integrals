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
| `exp_sub_one_sub_id_le_sq_mul_exp`, `norm_exp_sub_one_sub_id_le_sq` | `LocalLemma.lean` | $e^r-1-r\le\tfrac12r^2e^r$ for $r\ge0$; $\|e^w-1-w\|\le\tfrac12\|w\|^2e^{\|w\|}$ for complex $w$ | Lemma L (008) |
| `local_pointwise` | `LocalLemma.lean` | $e^{-\beta}\|e^w-1-q\|\le e^{-\beta}c_4+\tfrac12\|w\|^2e^{-(1-\kappa)\beta}$ when $\|w\|\le\kappa\beta$, $\|w-q\|\le c_4$ | Pointwise step of Lemma L (008) |
| `exp_mul_le_max_one_exp`, `hasDerivAt_ofReal_mul_const`, `hasDerivAt_cexp_ofReal_mul`, `norm_cexp_ofReal_mul_le`, `norm_cexp_ofReal_mul_sub_one_le`, `norm_cexp_sub_one_sub_id_le`, `norm_cexp_sub_one_sub_id_le_of_re_le` | `ExpRemainder.lean` | $\|e^w-1-w\|\le\tfrac12\|w\|^2\max(1,e^{\mathrm{Re}\,w})$ for complex $w$, and the form with $\mathrm{Re}\,w\le b$ | Second order remainder step of Lemma E-local (correspondence 049). Compiled in CI on 2026-11-24, pull request 60, axioms `propext`, `Classical.choice`, `Quot.sound`. The moment inequalities and the change of Gaussian measure of that lemma are not formalized. |

## Verified by interval arithmetic, not in Lean

| Item | Where | Remark |
| --- | --- | --- |
| Hypotheses of Theorem S at one parameter point, $\hbar=1/40,1/80,1/160$ | `experiments/2026-11-23-sparks/certify_polygon.py`, correspondence 010 | Arb ball arithmetic; program run by Sparks only |
| Hypotheses of Theorem S for three real saddles, symmetric ($c=0$) and asymmetric ($c=1/5$), $\hbar$ down to $1/320$ | `experiments/2026-11-23-sparks/certify_contour.py`, correspondence 016 | Arb ball arithmetic; program run by Sparks only |
| Two-variable selection of one complex saddle for a classically forbidden amplitude, $\hbar=1/20$ to $1/160$ | `experiments/2026-11-23-sparks/complex_selection/certify_2d.py`, correspondence 018 | Arb ball arithmetic; program run by Sparks only; depends on Lemma 1ₙ and Lemma Lₙ, unreviewed |
| Spectral quadrature budget and lattice enclosure at $N=8$ | Nova, `pa_interval.py`, correspondence 014 | mpmath interval backend; reproduced by Sparks |

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
| Propositions F and T, tilted-plane selection; separable certificate; homotopy with regulator | 019, 022, Nova 020 and 023 | Each reviewed by the other author |
| Convergence of the selection margin to a Sturm–Liouville eigenvalue; convergence of discrete complex saddles | Nova 034, 035 | Sparks (reviews on pull requests 39, 40) |
| Exterior and local certificates on the one-mode surface | Nova 039 to 046 | Sparks (047), partial: three derivations checked by hand |
| Lemma E-local, local error with N-independent constants | 049 | Under review by Nova (pull requests 57, 58); one step in Lean |

## Next Lean targets, in order

Priority follows the current line of work (correspondence 019 onward). Items 1 to 4 are short and
self-contained; item 5 is the first that needs Gaussian measures in Mathlib.

1. The matrix identity behind the tilted homotopy (Nova 023): $B_s=(1-s)D+sB_1+s(1-s)TDT$ and hence $B_s\ge sB_1$.
2. The endpoint trace identity and bound of Nova 045: $v_0+v_N-2\mu=\sum_j\frac{2j+1-N}{N}\Delta_j$ and $|v_0+v_N-2\mu|^2\le\tfrac23(1-N^{-2})E$.
3. The edge inequality for clipped slopes of Nova 039.
4. The damping inequalities of Nova 024: $t-\sin t\ge\min(t^3,t)/8$ and the convex minorant.
5. The two moment inequalities of Lemma E-local for a complex linear form of a Gaussian vector, $E|\delta|^4\le3(E|\delta|^2)^2$ and $E|\delta|^6\le15(E|\delta|^2)^3$.
6. The Gaussian moment integrals of Lemma L (its pointwise core is done).
7. The one-dimensional regulated Fresnel integral and its modulus bound.
8. Eigenvalues of the tridiagonal Hessian and the link to `Matrix.det`.
