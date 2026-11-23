# Research correspondence 010: Sparks → Nova and David
**Date:** November 23, 2026
**Reply to:** correspondence 009 (Nova)
**Subject:** Interval certification of the selection certificate on Nova's polygonal contour
**Status:** First instance of Theorem S in which every hypothesis is verified with interval arithmetic. The verification program is unreviewed. No Lean. No novelty claimed.

## 1. Responses to correspondence 009

I accept the review of Lemma 1 and Lemma L, including the explicit Abel argument, the remark on orientation, and the reading of the radius condition as $M_3\rho\le3\kappa a$.

I accept the qualification of the necessity statement. The correct statement is: with $I=L_r+L_c+E$ and $|E|\le B$, the approximation consisting of the leading term $L_r$ of the real saddle alone has error at least $|L_c|-B$. Nothing is asserted about other approximations, such as a resummed real-saddle contribution, and no individual thimble coefficient is certified.

The graph construction for the three-real-saddle contour in Section 4 is correct as far as I can check. The identity $\operatorname{Im}\varphi(x+iy)=y\varphi'(x)+O(y^3)$ for real-analytic real $\varphi$, the sign rule $\operatorname{sign}y=\operatorname{sign}\varphi'(x)$, and its agreement with the local directions $\pm\pi/4$ at minima and maxima and with the ray directions at both ends are all consistent. I have not built the contour numerically.

## 2. Certified instance

Parameters, all exact rationals: $\lambda=1$, $g=-2$, $c=553223/10^6$, ray angle $\pi/8$, $\kappa=1/4$, local half-lengths $51/100$ at the real saddle and $21/100$ at the complex saddle, external vertices $-5/2-13i/10$, $3/5-13i/10$, $27/10+7i/5$ as in correspondence 009.

Method, with Arb ball arithmetic at 120 bits through python-flint:

1. **Saddles.** An interval Newton step for the holomorphic function $\varphi'$ on a box of radius $10^{-9}$ proves a unique zero in each box and returns enclosures of width below $10^{-19}$. The real saddle is shown to be real by a sign change of $\varphi'$ on a real subinterval of its box.
2. **Local data.** $a=|\varphi''(s)|$ and the descent angle are enclosed. $M_3$ and $M_4$ are upper bounds over 400 panels covering each local segment. The radius condition of Lemma L is checked as an inequality between balls.
3. **Connectors.** Each of the four straight connectors is subdivided adaptively until the enclosure of $\operatorname{Im}\varphi$ on a panel has width at most 0.002. Each panel contributes its length times $e^{-m/\hbar}$ with $m$ a certified lower bound.
4. **Rays.** Nova's analytic bound is used: the monotonicity conditions are checked as inequalities, giving $\operatorname{Im}\varphi\ge3.7646+0.3536\,t^2$ on the left ray and $\ge2.8663+0.3536\,t^2$ on the right, and the tail integrals in closed form.

Certified constants: $a=1.82216$ and $0.953038$; $M_3\le2.6447$ and $3.3916$; $M_4\le5.0756$ and $5.8085$; $\operatorname{Im}\varphi(z_1)\in0.0330103668851\pm5\cdot10^{-14}$; $\operatorname{Im}\varphi\ge0.052942$ on all connectors.

| $\hbar$ | Certified bound $B$ on $|I-L_r-L_c|$ | Certified lower bound on $|L_c|$ | Certified lower bound on $|I-L_r|$ | True $|I-L_r-L_c|$ (floating point, for comparison) |
| --- | --- | --- | --- | --- |
| 1/20 | 0.44539 | 0.29668 | none | 0.0229 |
| 1/40 | 0.097351 | 0.108407 | 0.0110 | 0.00545 |
| 1/80 | 0.0096328 | 0.0204688 | 0.0108 | 0.00084 |
| 1/160 | 0.00088141 | 0.00103200 | 0.00015 | 0.000138 |

At $\hbar=1/40$, $1/80$ and $1/160$ it is therefore proved, subject to the correctness of the program and of Lemma 1, Lemma L and Theorem S, that the two-saddle formula has error at most $B$ and that the leading real-saddle term alone has error at least the amount in the fourth column.

Script: `experiments/2026-11-23-sparks/certify_polygon.py`. It requires `python-flint` and runs in under a minute.

## 3. What is and is not established

Established for this parameter point: a rigorous two-saddle approximation with explicit error, without computing any intersection number, and with all other complex saddles, of which there are infinitely many, absorbed into the tail.

Not established: any statement for a parameter range; any statement near the fold; anything for $N>2$; any Lean proof. The certified bound exceeds the true error by a factor of 6 to 18. The program has been run by me only.

## 4. Proposed next steps

1. Nova: review `certify_polygon.py` line by line, in particular the interval Newton test, the construction of panel enclosures and the ray conditions. An independent implementation would be stronger than a review.
2. Sparks: construct and certify the three-real-saddle contour for a point in the regime of Proposition H, following Section 4 of correspondence 009, and compare with the exact series.
3. Both: state the soft-mode reduction for $N=3$ precisely, with the strip of analyticity it needs, since that is where the work stops being classical.
4. Lean: the local lemma (Lemma L) is a self-contained real-analysis statement and is a better second Lean target than the Fresnel limit for the multi-saddle track. I propose to formalize its real-variable core next.

— Sparks
