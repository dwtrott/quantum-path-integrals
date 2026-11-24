# Research correspondence 018: Sparks → David and Nova
**Date:** November 23, 2026
**Subject:** Change of direction. Certified selection of a complex classical path in two variables: a classically forbidden amplitude
**Status:** Two lemmas with conventional proofs (unreviewed), one interval-certified instance (program unreviewed), floating-point exploration. No Lean. No novelty claim; the prior art named in Section 7 is cited from memory and is unaudited.

## 0. Why the direction changes

David's question is the original one: in what sense do only some paths contribute, and can that be proved. The work so far answers it in one variable and, in the single-saddle regime, repeats what Hargé proved. The part of the question that is open is the selection of complex classical paths in more than one variable. This note takes the first step there. It does not follow the numerical-analysis track of correspondences 014 to 017, which I leave to Nova.

## 1. Model: two kicks between wavepackets

Units $m=1$, time step $1$. A Gaussian wavepacket centred at $(a,p_0)$ propagates freely, is kicked by $e^{-igW(x)/\hbar}$ with $W(x)=e^{-x^2/2}$, propagates, is kicked again, propagates, and is projected on a Gaussian wavepacket centred at $(b,p_b)$. Both packets have width $\sqrt\hbar$. Integrating the two Gaussian end variables exactly leaves
$$I=\int_{\mathbb R^2}e^{i\Phi(x_1,x_2)/\hbar}dx_1dx_2,\qquad\Phi=\phi_{\rm in}(x_1)+\tfrac12(x_2-x_1)^2+\phi_{\rm out}(x_2)-g\big(W(x_1)+W(x_2)\big),$$
$$\phi_{\rm in}(x)=\tfrac12(x-a)^2+\tfrac{i-1}4\big(p_0-(x-a)\big)^2,\qquad\phi_{\rm out}(x)=\tfrac12(x-b)^2+\tfrac{i-1}4\big(p_b+(x-b)\big)^2 .$$
On $\mathbb R^2$, $\operatorname{Im}\Phi=\tfrac14[(x_1-a-p_0)^2+(x_2-b+p_b)^2]$, so the integral converges absolutely and needs no regulator. All coefficients are complex and, generically, so are all saddles.

Parameters: $a=-2$, $p_0=1$, $g=3$. Classically the packet is reflected: it reaches the first kick at $x=-1$ with momentum $1$, leaves it with momentum $-0.82$, and ends near $(b,p_b)=(-3.68,-1.86)$. The target studied here is the transmitted packet $(b,p_b)=(2,1)$, which is classically forbidden.

## 2. What the exploration shows

A Newton search in a box of $\mathbb C^2$ finds 123 saddles for the forbidden target. Sixty of them have $\operatorname{Im}\Phi<0$, down to $-54$, and would each contribute a term of size $e^{54/\hbar}$ if included. Sixty-three have $\operatorname{Im}\Phi\ge0$. The exact amplitude, computed by the trapezoid rule on the real plane, is reproduced by one of them alone,
$$z_s=(-0.284121-0.069441i,\;0.284121+0.069441i),\qquad\operatorname{Im}\Phi(z_s)=0.258118,$$
with relative error $3.1$, $1.5$ and $1.1$ percent at $\hbar=0.2,0.1,0.07$. The next saddles have $\operatorname{Im}\Phi\approx2.42$.

For the classically allowed target $(-3.68,-1.86)$ the leading saddle is nearly real, and a second saddle has $\operatorname{Im}\Phi=0.165>0$. Its term would be a correction of 10 to 40 percent at these values of $\hbar$, and the data show that it is absent: the error of the leading saddle alone is 3 to 8 percent and decreases linearly in $\hbar$. So a positive imaginary action does not imply a contribution. This is the selection problem in its plain form, and choosing saddles by the size of their action fails in both directions.

## 3. Two lemmas

**Lemma 1ₙ (deformation).** Let $\Phi=Q-g\sum_jW(z_j)$ with $Q$ a quadratic polynomial whose Hessian has positive definite imaginary part, and let $y:\mathbb R^n\to\mathbb R^n$ be bounded and Lipschitz. Then $\int_{\mathbb R^n}e^{i\Phi(x)/\hbar}dx=\int_{\mathbb R^n}e^{i\Phi(x+iy(x))/\hbar}\det(I+iDy(x))\,dx$.

*Proof.* The form $e^{i\Phi/\hbar}dz_1\wedge\dots\wedge dz_n$ is closed. Apply Stokes' theorem to the chain $\{x+ity(x):|x|\le r,\,0\le t\le1\}$. For bounded $y$, $|W(z_j)|\le e^{y_j^2/2}$ is bounded and $\operatorname{Im}Q(x+ity)\ge c|x|^2-C|x|-C'$ uniformly in $t$, so the lateral boundary term vanishes as $r\to\infty$. ∎

**Lemma Lₙ (local contribution).** Let $z_s$ be a nondegenerate saddle, $R=(-i\Phi''(z_s))^{-1/2}$ a symmetric root, so that $i\Phi(z_s+Ru)=i\Phi(z_s)-\tfrac12|u|^2+\hbar\,w(u)$. Let $M_3,M_4$ bound the third and fourth derivatives of $u\mapsto\Phi(z_s+Ru)$, as multilinear forms in the Euclidean norm, on the square $|u|_\infty\le\rho$, and assume $\sqrt n\,\rho\,M_3\le3\kappa$ with $0<\kappa<1$. Then
$$\int_{|u|_\infty\le\rho}e^{i\Phi(z_s+Ru)/\hbar}\det R\,du=(2\pi\hbar)^{n/2}\det R\;e^{i\Phi(z_s)/\hbar}(1+r),$$
$$|r|\le n\,\operatorname{erfc}\!\Big(\frac{\rho}{\sqrt{2\hbar}}\Big)+\hbar\Big[\frac{n(n+2)}{24}M_4+\frac{n(n+2)(n+4)}{72}\frac{M_3^2}{(1-\kappa)^{(n+6)/2}}\Big].$$

*Proof.* As for Lemma L in correspondence 008: the cubic term is odd and integrates to zero over the symmetric square; the pointwise inequality is `local_pointwise`, which is verified in Lean; the constants are the fourth and sixth moments of an $n$-dimensional Gaussian, $n(n+2)\hbar^2$ and $n(n+2)(n+4)\hbar^3$. For $n=1$ this is Lemma L. ∎

For the present $\Phi$ the higher derivatives come only from the kicks, and $M_k\le g\sum_j\sup|W^{(k)}(z_j)|\,\|r_j\|^k$ with $r_j$ the rows of $R$.

## 4. The surface

$$\Gamma=\{\,x_s+\operatorname{Re}(R)u+i\,(y_s+\operatorname{clip}_B(\operatorname{Im}(R)u)):u\in\mathbb R^2\,\},$$
where $z_s=x_s+iy_s$ and $\operatorname{clip}_B$ truncates each component to $[-B,B]$. Near the saddle $\Gamma$ is the exact tangent plane of the thimble, so Lemma Lₙ applies on the square. Away from it $\Gamma$ stays in a strip of half-width about $B$, where the Gaussian kick is bounded. It is piecewise affine, so its Jacobian is piecewise constant. The certificate is then the statement of Theorem S: the integral equals the local term plus the integral of $|e^{i\Phi/\hbar}|$ over the rest of $\Gamma$.

## 5. Certified result

Program: `experiments/2026-11-23-sparks/complex_selection/certify_2d.py`. Arb ball arithmetic, 100 bits. Parameters exact rationals; $B=3/5$, $\rho=9/20$, $\kappa=1/4$.

1. The saddle is certified by a Krawczyk test in $\mathbb C^2$ on a box of radius $10^{-8}$, with enclosure width below $10^{-13}$.
2. $R$ is computed from the closed form of the $2\times2$ square root; $\det\operatorname{Re}R>0$ is certified, which fixes the orientation and the sign of the leading term.
3. $M_3\le1.0846$, $M_4\le1.1611$, and $\sqrt2\rho M_3=0.690<0.75$.
4. The surface outside the square is covered by 2.6 million adaptive panels on $|u|_\infty\le48$; on every panel $\operatorname{Im}\Phi-\operatorname{Im}\Phi(z_s)\ge0.0973$.
5. Beyond that, an explicit bound $\operatorname{Im}\Phi\ge\sum_js_j^2/4-1.674\sum_j|s_j|-11.75$, with $s_j$ the distance from the envelope centre, gives a closed-form Gaussian tail.

| $\hbar$ | Modulus of the leading term | Certified bound on $|I-L_0|/|L_0|$ | of which local | of which surface tail | True relative error (floating point) |
| --- | --- | --- | --- | --- | --- |
| 1/20 | $4.43\cdot10^{-4}$ | 0.322 | 0.232 | 0.091 | 0.0076 |
| 1/40 | $1.27\cdot10^{-6}$ | 0.0901 | 0.0805 | 0.0096 | 0.0038 |
| 1/80 | $2.08\cdot10^{-11}$ | 0.0361 | 0.0359 | 0.00013 | 0.0019 |
| 1/160 | $1.12\cdot10^{-20}$ | 0.0179 | 0.0179 | $3\cdot10^{-8}$ | not computable in double precision |

Subject to Lemma 1ₙ, Lemma Lₙ and the correctness of the program, this proves that a classically forbidden amplitude of size $10^{-20}$ is given, to within 1.8 percent, by the contribution of one complex classical path. All other saddles in $\mathbb C^2$, whatever their number and their intersection numbers, are inside the bound. In particular none of the saddles with negative imaginary action contributes, and none of the others contributes more than the stated fraction.

The run takes about one minute.

## 6. What is not done

1. **The exclusion case.** For the classically allowed target I tried the same surface around the nearly real saddle, with the aim of proving that the competitor with $\operatorname{Im}\Phi=0.165$ is absent. The design fails: the tangent plane there tilts by about 45 degrees and passes through regions where $\operatorname{Im}\Phi$ is negative, and $M_3\approx8.5$ forces $\rho\approx0.06$. A different surface is needed, probably one that follows the real plane more closely. This is the next target, because a certified exclusion is the stronger form of selection.
2. **One parameter point.** No range of targets or of $\hbar$. No Stokes transition has been located or crossed.
3. **Two variables only.** The panel count grows exponentially with dimension, so this method by itself stops at a few variables. The reduction of correspondences 011 to 013 is what could extend it.
4. **No rule.** The result certifies a selection in an instance. It does not say why this saddle and not another. A rule that predicts the contributing set for a family, and a proof of it, is the result that would matter, and I have no candidate yet.
5. The lemmas are unreviewed and nothing here is in Lean beyond the pointwise inequality.

## 7. Prior art to audit

From memory and unverified: Shudo and Ikeda on complex classical trajectories and chaotic tunneling in kicked systems, where the selection of contributing complex trajectories by Stokes geometry is the central difficulty; Huber, Heller and Littlejohn on generalized Gaussian wavepacket dynamics; work on Stokes phenomena for multidimensional saddle point integrals (Howls, and Delabaere and Howls); the generalized thimble and holomorphic flow methods in lattice field theory; Feldbrugge and Turok on real-time path integrals by Picard–Lefschetz theory. I expect certified error bounds for a multi-step kicked system to be absent from this literature, and I have not checked.

## 8. Requests

To Nova: review Lemma 1ₙ and Lemma Lₙ, in particular the moment constants and the claim that the lateral boundary term vanishes; and review or reimplement `certify_2d.py`. To David: Sections 2 and 3 are the parts where a reader with your background adds most.

— Sparks
