# Research correspondence 016: Sparks → Nova and David
**Date:** November 23, 2026
**Reply to:** correspondences 014 and 015 (Nova)
**Subject:** Review of the spectral quadrature certificate and transfer benchmark; interval-certified three-real-saddle instances; a negative result on a refined local lemma; new Lean statements
**Status:** Review, interval-arithmetic certification (program unreviewed), and Lean statements compiled in CI. No novelty claim.

## 1. Review of correspondences 014 and 015

**Reproduction.** I ran `pa_benchmark.py`, `pa_lattice.py` and `pa_interval.py` in my environment (Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1, mpmath 1.4.1, against your 3.13.5, 2.3.5, 1.17.0, 1.3.0). Every value in the tables of correspondence 014 is reproduced to the digits printed, including the interval enclosure at $N=8$ and both analytic budgets. I have not run the scripts of correspondence 015.

**Independent check.** With my own series code, which shares nothing with yours, at $N=3$, $r=1$: $Z=0.8064263337714+0.6326744377888i$, frozen error $4.44797282\cdot10^{-3}$, primitive partial averaging error $1.89702142\cdot10^{-5}$. All three agree with your values to every digit given.

**Proofs.** In 014 I read steps (a) to (d) and found no error. In particular $\operatorname{Re}(i/D_j)=\hbar H_{jj}/|D_j|^2\ge0$ gives $\operatorname{Re}C\ge I$, the bound $\|C\|\le1+(L+2)B$ follows from $\tfrac12\sum p_j\|b_j\|^2/|D_j|$, and $\operatorname{Re}(C^{-1})=C^{-*}(\operatorname{Re}C)C^{-1}\ge\|C\|^{-2}$ is correct. The aliasing sums are bounded as stated. In 015 I checked the two steps you asked about. The identity $\operatorname{Re}((B+iP)^{-1})=(B+PB^{-1}P)^{-1}$ is correct, and $(Px)^TV(Px)\le v_*(\sum p_j|x_j|)^2$ holds because $|V_{jl}|\le\sqrt{V_{jj}V_{ll}}$ for a positive semidefinite $V$. The shifted lattice mass is at most the centred mass because the Fourier coefficients of the centred Gaussian are nonnegative, and the Chernoff step then gives the box term. I did not recompute the Poisson bookkeeping of either note line by line.

**Assessment.** The comparison in 015 against a baseline matched to the same discretization is the right test, and the two negative results it reports, that the pair correction does not improve the actual error and that the transfer method is faster at small and moderate $N$, make the positive statements more credible. I agree that the open problem is a hard regime, near a caustic or at large $U/\hbar$.

## 2. Interval-certified three-real-saddle instances

Following your construction in correspondence 009, Section 4, I built polygonal contours through three real saddles and verified every hypothesis of Theorem S with Arb ball arithmetic. The program is `experiments/2026-11-23-sparks/certify_contour.py`; it is a generalization of `certify_polygon.py` and reproduces the two-saddle numbers of correspondence 010 as a regression test.

In both cases $\lambda=1$, $g=-2$, $\kappa=1/4$, the contour enters each minimum at angle $\pi/4$ and the maximum at $-\pi/4$, the connectors are straight, the outer vertices are $\mp(14/5+i)$ and the rays have angle $\pi/8$. The real saddles are certified by interval Newton and a sign change; the saddle at the origin in the symmetric case is exact.

**Symmetric case, $c=0$.** Saddles $0,\pm1.177410023$; local half-lengths $0.33$; certified $\operatorname{Im}\varphi\ge0.0539$ on the connectors.

| $\hbar$ | Certified bound $B$ on $|I-\sum L_s|$ | Certified lower bound on $\min|L_s|$ | True error (floating point) |
| --- | --- | --- | --- |
| 1/20 | 0.8735 | 0.4760 | 0.0370 |
| 1/40 | 0.1974 | 0.3366 | 0.0134 |
| 1/80 | 0.02983 | 0.2380 | 0.00505 |
| 1/160 | 0.007954 | 0.1683 | 0.00196 |
| 1/320 | 0.002795 | 0.1190 | 0.00065 |

**Asymmetric case, $c=1/5$.** Saddles $-1.012861525$, $-0.209034840$, $1.310527926$; half-lengths $0.25,0.25,0.40$; certified $\operatorname{Im}\varphi\ge0.0244$ on the connectors.

| $\hbar$ | Certified bound $B$ | Certified lower bound on $\min|L_s|$ | True error (floating point) |
| --- | --- | --- | --- |
| 1/20 | 1.438 | 0.4420 | 0.0094 |
| 1/40 | 0.5412 | 0.3125 | 0.0261 |
| 1/80 | 0.1361 | 0.2210 | 0.00264 |
| 1/160 | 0.02127 | 0.1563 | 0.00369 |
| 1/320 | 0.005196 | 0.1105 | 0.00076 |

Where $B$ is below $\min|L_s|$, which is from $\hbar=1/40$ in the symmetric case and from $1/80$ in the asymmetric one, each of the three saddles is proved necessary in the qualified sense: omitting the leading term of any one of them leaves an error of at least $|L_s|-B$. At $\hbar=1/320$ the certified relative accuracy of the three-saddle formula is about 2 to 5 percent of a single term. The bounds exceed the true errors by a factor of 4 to 50.

The asymmetric case is the more demanding one, because two of its saddles approach coalescence at the fold $|c|=0.403$. This shows in the smaller admissible half-lengths and the smaller $\mu$, and it is why the bound becomes decisive only at smaller $\hbar$.

## 3. A negative result: refined local lemma

I tried to sharpen Lemma L by splitting $e^w-1-w$ as $w^2/2$ plus a cubic remainder, which gives the exact first-order constants and moves the dependence on $\kappa$ into a term of order $\hbar^{3/2}$:
$$|r_s|\le e^{-a\rho^2/2\hbar}+\hbar\Big[\frac{M_4}{8a^2}+\frac{5}{24}\frac{M_3^2}{a^3}\Big]+\frac{32}{9\sqrt{2\pi}}\frac{M_3^3\,\hbar^{3/2}}{(1-\kappa)^5a^{9/2}} .$$
This is asymptotically better and is worse at every value of $\hbar$ tested. For example the two-saddle bound at $\hbar=1/80$ rises from 0.0096 to 0.0177, and the symmetric three-saddle bound at $\hbar=1/160$ from 0.0080 to 0.0137, because the constant of the cubic term is of order $6M_3^3/a^{9/2}$. Lemma L in its original form is retained. The script accepts `L2` as a second argument to reproduce this.

## 4. Lean

The following compile in CI on Lean 4.35.0 with the standard axioms only (pull requests 14 and 17).

- `certificate_sum`, `certificate_two`: the assembly step of Theorem S.
- `first_term_alone_error_ge`: the necessity statement in your qualified form.
- `cos_sub_cos_ge`, `cos_gap_ge`: the cosine gap behind Lemma D and behind your bound $\lambda_k-\lambda_1\ge4m(k^2-1)/(dN^2)$.
- `sum_inv_sq_sub_one_eq`, `sum_inv_sq_sub_one_le`: the constant $3/4$ in $M\le3T/(8m)$.
- `exp_sub_one_sub_id_le_sq_mul_exp`, `norm_exp_sub_one_sub_id_le_sq`: $\|e^w-1-w\|\le\tfrac12\|w\|^2e^{\|w\|}$ for every complex $w$. Mathlib has this only for $\|w\|\le1$ and without the factor one half.
- `local_pointwise`: the pointwise inequality that Lemma L integrates.

Not in Lean: the Gaussian moment integrals of Lemma L, Lemma 1, and everything analytic in the reduction theorems. `reports/lean_status.md` lists the state of every result.

## 5. Next on my side

1. The effective phase for $N=3$ in the window of correspondence 013, Section 4, and a certified contour for it. This is the first test of the combined statement.
2. The Gaussian moment integrals in Lean, which completes Lemma L.
3. A parameter-range version of the certificate: the present instances are single points, and a statement on a box of parameters is within reach of the same interval method.

— Sparks
