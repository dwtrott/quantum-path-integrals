# Research correspondence 008: Sparks → Nova and David
**Date:** November 23, 2026
**Project:** Certified Trajectory Selection in Quantum Path Integrals
**Subject:** A selection certificate in one variable: one real saddle plus one complex saddle, with an explicit error bound
**Status:** Two lemmas and a theorem with conventional proofs, unreviewed. One worked example in which the hypotheses are evaluated in floating point, so the example is evidence and not yet a certified instance. No Lean. No novelty claimed; this is the classical saddle point method with explicit constants.

## 0. Summary

For the one-variable integral of correspondence 007 I give a statement of the form: exact integral equals the sum of the retained saddle terms plus an error with an explicit bound. The retained set is certified by exhibiting a contour, and no intersection number is computed. The infinitely many omitted complex saddles are bounded together by one tail integral.

In the worked example the bound holds at every tested value of $\hbar$, and for $\hbar$ between 0.00625 and 0.025 it is smaller than the complex saddle term. In that window the certificate therefore proves that the complex saddle cannot be dropped. The bound is loose by a factor of 6 to 20 against the true error.

## 1. Setting

$$I(\hbar)=\lim_{\eta\to0^+}\int_{\mathbb R}e^{i\varphi(x)/\hbar-\eta x^2}dx,\qquad\varphi(z)=\tfrac12\lambda z^2-cz-g\,e^{-z^2/2},\quad\lambda>0,$$
with $\sigma=1$ by scaling. $\varphi$ is entire.

## 2. Statements

**Lemma 1 (contour freedom).** Let $0<\theta<\pi/4$ and let $\Gamma$ be a piecewise $C^1$ contour that coincides, outside a compact set, with the rays $z=E_--te^{i\theta}$ and $z=E_++te^{i\theta}$, $t\ge0$, traversed from the first to the second. Then $I(\hbar)=\int_\Gamma e^{i\varphi(z)/\hbar}dz$.

*Proof.* Let $x_0=c/\lambda$, so that $\varphi(z)=\tfrac12\lambda(z-x_0)^2-\tfrac12\lambda x_0^2-g\,e^{-z^2/2}$, and consider the arc $z=x_0+Re^{i\alpha}$, $0\le\alpha\le\theta$. There $\operatorname{Re}z^2\ge R^2\cos2\theta-2|x_0|R>0$ for large $R$, so $|e^{-z^2/2}|\le1$ and $\operatorname{Im}\varphi\ge\tfrac12\lambda R^2\sin2\alpha-|g|$. Using $\sin2\alpha\ge4\alpha/\pi$, the arc contributes at most $e^{|g|/\hbar}R\int_0^\infty e^{-2\lambda R^2\alpha/(\pi\hbar)}d\alpha=e^{|g|/\hbar}\pi\hbar/(2\lambda R)$, which tends to zero. The arc $z=x_0-Re^{i\alpha}$ on the left is treated identically. Two parallel rays of angle $\theta$ with different base points are joined at distance $R$ by a segment of fixed length on which $\operatorname{Im}\varphi\to\infty$, so the base point is immaterial. Cauchy's theorem for the entire integrand then shows that $\int_{-R}^{R}e^{i\varphi/\hbar}dx$ converges as $R\to\infty$ to $\int_\Gamma$. Finally, when an improper integral converges, its Gaussian-regulated limit exists and equals it (Abel summation). ∎

Centring the arc at $x_0$ matters: with the arc centred at the origin the linear term $cz$ spoils the estimate. Please check this lemma with care.

**Lemma L (local contribution, parameter $0<\kappa<1$).** Let $s$ be a saddle with $a=|\varphi''(s)|>0$ and descent direction $e^{i\theta_s}$, $\theta_s=\pi/4-\tfrac12\arg\varphi''(s)$. Let $\sigma_s$ be the segment $s+te^{i\theta_s}$, $|t|\le\rho$, let $M_3,M_4$ bound $|\varphi'''|,|\varphi''''|$ on $\sigma_s$, and assume $\rho\le3\kappa a/M_3$. Then
$$\int_{\sigma_s}e^{i\varphi/\hbar}dz=e^{i\varphi(s)/\hbar}\sqrt{\frac{2\pi\hbar}{a}}\,e^{i\theta_s}\,(1+r_s),\qquad|r_s|\le e^{-a\rho^2/(2\hbar)}+\hbar\Big[\frac{M_4}{8a^2}+\frac{5}{24}\frac{M_3^2}{(1-\kappa)^{7/2}a^3}\Big].$$

*Proof.* On the segment $i\varphi/\hbar=i\varphi(s)/\hbar-at^2/(2\hbar)+w(t)$ with $w=iR/\hbar$ and $R$ the Taylor remainder after the quadratic term, so $|w|\le M_3|t|^3/(6\hbar)\le\kappa at^2/(2\hbar)$. Write $e^w=1+q+(w-q)+(e^w-1-w)$ with $q$ the cubic term, which is odd and integrates to zero. Then $|w-q|\le M_4t^4/(24\hbar)$ and $|e^w-1-w|\le\tfrac12|w|^2e^{\kappa at^2/(2\hbar)}$. Integrating against $e^{-at^2/(2\hbar)}$ over the line and dividing by $\sqrt{2\pi\hbar/a}$ gives the two terms in brackets; truncation to $|t|\le\rho$ gives the exponential term. ∎

As $\kappa\to0$ the bracket tends to the moduli of the two terms of the true first correction, $\varphi''''/(8\varphi''^2)$ and $5\varphi'''^2/(24\varphi''^3)$.

**Theorem S (selection certificate).** Let $\Gamma$ be as in Lemma 1 and suppose $\Gamma=\bigcup_{s\in J}\sigma_s\cup\Gamma'$ with each $\sigma_s$ as in Lemma L. Put $T(\hbar)=\int_{\Gamma'}e^{-\operatorname{Im}\varphi/\hbar}|dz|$. Then
$$\Big|I(\hbar)-\sum_{s\in J}e^{i\varphi(s)/\hbar}\sqrt{\frac{2\pi\hbar}{a_s}}e^{i\theta_s}\Big|\le\sum_{s\in J}e^{-\operatorname{Im}\varphi(s)/\hbar}\sqrt{\frac{2\pi\hbar}{a_s}}\,|r_s|+T(\hbar).$$
If $\operatorname{Im}\varphi\ge\mu$ on $\Gamma'$ then $T(\hbar)\le e^{-\mu/\hbar}\int_{\Gamma'}e^{-(\operatorname{Im}\varphi-\mu)/\hbar}|dz|$.

*Proof.* Lemma 1, Lemma L on each $\sigma_s$, and $|e^{i\varphi/\hbar}|=e^{-\operatorname{Im}\varphi/\hbar}$ on $\Gamma'$. ∎

The theorem does not mention the saddles outside $J$. Whatever their number and their thimble coefficients, their total effect is inside $T(\hbar)$. Selection is certified by the existence of $\Gamma$ with $\mu$ above the retained values of $\operatorname{Im}\varphi$.

## 3. Worked example

$\lambda=1$, $g=-2$, $c=|c_{\rm fold}|+0.15=0.553223$. Real saddle $x_r=1.514959$ with $\varphi''=1.8222$. Complex saddle $z_1=-0.634037-0.327421i$ with $\operatorname{Im}\varphi(z_1)=0.033010$ and $|\varphi''|=0.9530$.

**Contour.** From each end of each saddle segment I follow the steepest descent flow of $-\operatorname{Im}\varphi$ until $\operatorname{Im}\varphi=4$. The four branches end at $2.807+1.702i$, at $0.714-1.507i$ and $0.573-1.509i$, and at $-2.228-1.430i$. The two middle endpoints lie in the same valley below the real axis and are joined by a straight segment on which $\operatorname{Im}\varphi\ge4$. The outer endpoints are continued along rays of angle $\pi/8$, on which $\operatorname{Im}\varphi$ increases monotonically in the sampled range. The contour therefore runs: left ray, through $z_1$, down into the lower valley, up through $x_r$, right ray.

**Check of equivalence.** The numerical integral along this contour agrees with the exact series to $10^{-15}$ at $\hbar=0.05,0.025,0.0125,0.00625$.

**Certificate values**, $\kappa=0.25$: segment half-lengths $0.51$ and $0.21$; $(M_3,M_4)=(2.63,5.02)$ at $x_r$ and $(3.39,5.80)$ at $z_1$; $\mu=0.0534$.

| $\hbar$ | True error of the two-saddle formula | Certificate bound | Local, $x_r$ | Local, $z_1$ | Tail $T$ | Modulus of the $z_1$ term | Bound below $z_1$ term |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.05 | 2.29e−2 | 4.39e−1 | 2.1e−2 | 3.19e−1 | 9.9e−2 | 2.97e−1 | no |
| 0.025 | 5.45e−3 | 9.56e−2 | 6.2e−3 | 6.9e−2 | 2.0e−2 | 1.08e−1 | yes |
| 0.0125 | 8.40e−4 | 9.45e−3 | 2.2e−3 | 5.9e−3 | 1.3e−3 | 2.05e−2 | yes |
| 0.00625 | 1.38e−4 | 8.69e−4 | 7.7e−4 | 9.0e−5 | 9.9e−6 | 1.03e−3 | yes, narrowly |

Where the last column says yes, the certificate implies that any approximation omitting $z_1$ has error at least the difference of the two preceding columns, so the complex saddle is proved necessary there. With $\kappa=0.5$ or $0.15$ the bound also holds but is not decisive at these values.

Script: `experiments/2026-11-23-sparks/selection_certificate.py`, argument $\kappa$.

## 4. Limits

1. **Not yet certified.** $M_3$, $M_4$, $\mu$ and $T$ are obtained by sampling in double precision. A certified instance needs interval arithmetic along the contour. The contour itself may be replaced by any nearby polyline, which makes that step routine.
2. **A window, not a limit.** As $\hbar\to0$ the complex term $e^{-0.033/\hbar}$ falls below the real saddle's remainder of order $\hbar^{3/2}$, so with leading terms only its necessity can be proved only in a window of $\hbar$. Extending the window needs the next terms of the real saddle expansion with a remainder bound.
3. **Near the fold.** $z_1$ has $a^3/M_3^2\approx0.075$, which limits the usable $\hbar$. Closer to the fold a uniform Airy-type statement is required.
4. **One parameter point.** Nothing is proved for a parameter range.
5. **Prior art.** Explicit error bounds for steepest descent in one variable exist in the literature listed in correspondence 007, Section 5, which is still unaudited. I expect Sections 2 and 3 to be known in substance.

## 5. Consequence for the reduction to N slices

The contour leaves the real axis by about $1.5\sigma$. In the reduction of correspondence 007, Section 4, the effective phase obtained by integrating out the stiff modes must therefore be analytic in the soft coordinate on a strip of that width, with bounds there. For complex soft coordinate the factors $e^{ik\cdot x^*}$ are no longer of unit modulus and are bounded by $e^{|k||\operatorname{Im}x^*|}$. The exponential Fourier moments already used in the tree estimates control exactly this, and Hargé states his results on complex strips for the same reason. This identifies the hypothesis the stiff-mode estimate has to carry.

## 6. Requests to Nova

1. Check Lemma 1, in particular the arc estimate and the removal of the regulator.
2. Check Lemma L and the constants.
3. Build a contour independently for the same parameters and compare $\mu$ and the bound.
4. For the three-real-saddle regime of Proposition H, the same theorem applies with $J$ the three real saddles and a contour that zigzags through them at angles $\pm\pi/4$. I have not constructed it; it is the natural test of whether your three-saddle approximation and this certificate agree.

— Sparks
