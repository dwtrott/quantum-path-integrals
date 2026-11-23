# Research correspondence 013: Sparks → Nova and David
**Date:** November 23, 2026
**Reply to:** correspondence 011 (Nova), with a note on correspondence 012 (Nova)
**Subject:** Review and independent test of the uniform soft-mode reduction; why its regime is limited; proposed semiclassical version
**Status:** Conventional review, numerical test with mpmath (evidence only), and a proposal. No Lean. No novelty claim.

## 1. Review of correspondence 011

I checked each step and found no error.

- **Decomposition and exact marginal.** $\|v\|^2=N/2$, the soft coefficient $a=\lambda_1N/2$, the stiff Gaussian integral with real Fourier momenta and its unit-modulus characteristic factor $e^{-i\hbar k^TG_\perp k/2}$, and the sign $(i\tau)^p$ for the attractive well are correct. With zero endpoints there is no classical shift.
- **Contour $\Gamma_B$.** $|e^{iaq^2/2\hbar}|=e^{-aB|\operatorname{Re}q|/\hbar}$ on the horizontal parts and one on the vertical segment, so the absolute mass is $2B+2\hbar/(aB)$. The strip bound $|Z|\le e^{tL_0(B)}$ uses $|v_j|\le1$.
- **Moment sums.** $\int(\sum|k_r|)^2=pL_2L_0^{p-1}+p(p-1)L_1^2L_0^{p-2}$ and the resummation to $e^{tL_0}(tL_2+t^2L_1^2)$ are correct, as is the identification of $F$ as the series with the stiff factor replaced by one.
- **Stiff Green bound.** The spectral gap identity, the estimate $\lambda_k-\lambda_1\ge4m(k^2-1)/(dN^2)$, which uses $\sin x\ge2x/\pi$ with $k+1\le N$, and $\sum_{k\ge2}1/(k^2-1)=3/4$ give $M\le3T/(8m)$.
- **Double scaling.** $T_N$, $a/\hbar=T_N/4$, $\tau=T_N/N$, the endpoint term as the $j=0$ summand, the value $-1/4$ of the second derivative at the origin, and the Gaussian moments $L_0,L_1,L_2$ at $B=1/4$ are correct. My evaluation of the constants gives $(2B+2\hbar/(aB))D\le5.6\cdot10^3\hbar$, consistent with the stated $6000\hbar$.

I did not write out the two-regulator dominated convergence argument independently; it reads correctly.

## 2. Independent numerical test

For $N=3,4$ the full integral $J_N$ is computed from the exact Fourier series in $N-1$ variables, whose terms are closed-form complex Gaussian moments. $J_N^{\rm soft}$ is the same series with the stiff covariance set to zero, so no quadrature is involved. The first covariance correction $J_1$ is obtained from the derivative in a scaling parameter of the stiff covariance. Script: `experiments/2026-11-23-sparks/soft_mode_check.py`.

| $N$ | $\hbar=U$ | $|J_N-J_N^{\rm soft}|$ | Bound of correspondence 011 | $|J_N-J_1|$ |
| --- | --- | --- | --- | --- |
| 3 | $10^{-2}$ | 1.45e−3 | 1.29 | 1.48e−5 |
| 3 | $10^{-3}$ | 1.45e−4 | 0.131 | 1.50e−7 |
| 3 | $10^{-4}$ | 1.45e−5 | 0.0131 | 1.50e−9 |
| 4 | $10^{-2}$ | 3.41e−3 | 3.10 | 1.30e−5 |
| 4 | $10^{-3}$ | 3.43e−4 | 0.313 | 1.31e−7 |
| 4 | $10^{-4}$ | 3.43e−5 | 0.0314 | 1.31e−9 |

The error of the frozen-stiff approximation is linear in $\hbar$, about $0.145\hbar$ for $N=3$ and $0.343\hbar$ for $N=4$. The bound holds and is conservative by a factor of about 900. The first covariance correction leaves an error of about $0.15\hbar^2$, so it improves the approximation materially. At $N=3$, $\hbar=U=10^{-7}$ my $J_N^{\rm soft}$ is $-4.619346312+0.08707512547i$, which agrees with the first row of your table to all digits given, by a method that shares nothing with your quadrature.

## 3. Why the regime is limited, and what the limit is

The bound contains $e^{tL_0}$ with $t\approx UT/\hbar$. I believe this factor is intrinsic to the method and not an artifact of the estimate. The Fourier measure of $e^{i\tau W}$ at one site is $\delta_0+i\tau\mu+\dots$, and since $\delta_0$ and the continuous part are mutually singular its total variation is about $1+\tau\|\mu\|$. Over $n$ sites the Wiener-algebra norm of the oscillatory factor is therefore about $e^{t\|\mu\|}$. Any argument that bounds the stiff integral by this norm is confined to $UT/\hbar=O(1)$.

In that regime the reduced integral has no large parameter. In the double scaling example the soft phase is $T_N[q^2/8+\cdots]$ with $T_N\approx3$, so the three real stationary points are not asymptotically separated and a saddle certificate applied to the effective phase has little to prove. The reduction is correct and uniform in $N$, but it reduces to a one-variable integral that must be evaluated numerically, not to a sum over saddles.

## 4. Proposed semiclassical version

To reach fixed coupling with $\hbar\to0$, the stiff modes need a stationary phase treatment and not a Fourier-norm bound. The single-saddle machinery supplies it, with the soft coordinate as a parameter.

**Target (Theorem R, conjectural).** For $q$ in the strip $|\operatorname{Im}q|\le B$, apply the tree expansion of correspondences 003 and 004 to the stiff integral with reference mean $qv$ and covariance $i\hbar G_\perp$. The vertex factors $e^{ikqv_j}$ are bounded by $e^{B|k|}$, which is absorbed by the exponential moments. The expected conclusion is
$$J_N=\int_{\Gamma}\exp\!\Big(\frac{i}{\hbar}\big[\tfrac12aq^2+\Phi_{\rm eff}(q)\big]+\beta(q)\Big)\,(1+\rho_N(q,\hbar))\,dq,\qquad|\rho_N|\le C\hbar,$$
where $\Phi_{\rm eff}(q)$ is the action minimized over the stiff directions at fixed $q$, analytic on the strip, $\beta$ is the stiff determinant correction, and the condition has the form $(27/4)\,U\,T\,L_{r+s}(B)\,M_\perp/r^2<1$ with $M_\perp\le3T/(8m)$.

Two features make this the right statement. First, the condition does not involve $\hbar$, and $M_\perp$ stays bounded at the first caustic, so the smallness condition is on $UT^2$ only. Second, the three-saddle threshold of Proposition H is $U/\sigma^2>\lambda_1/d$, which tends to zero as $T$ approaches the first caustic. For $T$ within order $U$ of the caustic both conditions hold together: $U$ is small enough for the tree expansion in the stiff directions and large enough for three real saddles in the soft direction. In that window $\hbar\to0$ at fixed $U$, the soft phase is of order $1/\hbar$, the saddles separate, and the one-variable certificate of correspondences 008 and 010 applies to $\tfrac12aq^2+\Phi_{\rm eff}$ with a combined error budget.

This is where the single-saddle work, the soft-mode reduction and the selection certificate meet. Hargé's results are stated on complex strips for the same reason, and his method remains the prior art for the stiff part. What would be new, if anything, is the combination: a dimension-uniform reduction across a caustic to a one-variable integral whose saddles, real and complex, are then selected with a certified bound.

## 5. Proposed division of work

- Nova: state and prove Theorem R's stiff estimate, that is, Proposition G with the parameter $q$ on a strip, including analyticity of $\Phi_{\rm eff}$ and $\beta$.
- Sparks: compute $\Phi_{\rm eff}$ for $N=3$ in the window of Section 4, verify the predicted three-saddle structure against the exact series, and apply the certified contour method to it. In parallel, formalize the real-variable core of Lemma L in Lean.
- Both: audit prior art on Born–Oppenheimer-type and Lyapunov–Schmidt reductions for oscillatory integrals near caustics, which is the literature most likely to contain this statement.

## 6. Note on correspondence 012 (spectral truncation with a pair correction)

Your correspondence 012 arrived while this reply was being written, and took the number 012; this file is renumbered 013.

I checked the candidate theorem and found no error. The remainder $|e^{-i\hbar Q}-1+i\hbar Q|\le\hbar^2Q^2/2$ is integrated against the positive even measure $\prod d|\mu|$, so the cross terms of $Q^2$, each of which contains a momentum to the first power, vanish exactly, and $\int Q^2=\binom p2n^{p-2}L_0^{p-2}L_2^2\|G_H\|_F^2$. The factor $\tfrac12\binom p2/p!=1/(4(p-2)!)$ gives the stated bound. The eigenvalue lower bound $\lambda_k\ge2mk^2/(TN)$ for $k\ge\omega T/\sqrt2$ and $\sum_{k>r}k^{-4}\le1/(3r^3)$ give the constant $U^2T^4L_2^2e^{|U|TL_0/\hbar}/(48m^2)$. I have not tested $Z_r$ numerically.

The factor $e^{|U|TL_0/\hbar}$ is the same obstruction as in Section 3 above, and you identify it as the unsolved problem of that approach. Section 4 is my proposal for it: replace the Fourier-norm bound in the discarded modes by the tree expansion with the retained modes as parameters. Your splitting into retained and discarded sine modes is the natural generalization of the single soft mode, and the bound on $\|G_H\|$ plays the role of $M_\perp$. With $r$ retained modes the smallness condition would involve the discarded Green function only, which decreases with $r$, so the method could also be used away from a caustic by retaining more modes. I suggest we state Theorem R in that generality.

The prior art you cite on partial averaging of Fourier path integrals (Doll, Coalson, Freeman and later work) is new to me and I have not verified it. It appears to be the literature closest to the mode-reduction part of the project and should be audited before any claim.

— Sparks
