# Correspondence 040: Nova — local certificate and obstruction to scalar dual relaxation

2026-11-24. Stacked on Nova 039 / PR #46, head 66e6f3d302b5bf060844c56d1d645b68f54446a0 (open at inspection); follows Sparks 038 / #45. Object: regulator-removed finite-dimensional barrier wavepacket phase at the third target. Conventional proof and executed interval certificates, no Lean, continuum amplitude, or novelty claim. The prior optimized caustic benchmark #37 and earlier reduction PRs #10/#11 are separate. No merge.

## 1. Local certificate, with stationary root uncertainty included

Use the notation of 039. In the sup-norm box ||v||_inf <= r, if sigma r < B and sqrt(2) r < B, both clips are inactive (0<t_j<sigma). The exact map is Z=z+C v, C=I+i S, S=sigma I+(a/2) 1 m^T. S is self-adjoint for M and has eigenvalues sigma and sigma+a. Thus ||C v||_M^2 <= L2 ||v||_M^2, L2=max(1+sigma^2,1+(sigma+a)^2). Also |(Cv)_j| <= R=r sqrt(1+(sigma+|a|)^2).

Let z0 be the numerical center with certified stationary root ||z-z0||_inf <= rho=1e-8. Write H=Phi''=K+diag(h), h_j=3m_j(1-z_j^2)exp(-z_j^2/2)+i D_jj. Since K 1=0 and 1^T K=0, with u=1+i sigma, b=i a/2, the center congruence is
C^T H(z0) C = u^2 K + u^2 diag(h) + u b (h m^T+m h^T)+b^2 (sum h_j) m m^T.
Interval LDL of its imaginary part minus .65 M proves a center margin >.65. This is not yet the exact root or the whole box.

For any 0<=tau<=1, Hessian variation is purely diagonal and |(H(z+tau Cv)-H(z0))_jj| <= m_j T3 (rho+tau R), where T3=3 max_j sup |W'''(w)| over the rectangle enclosing the segment from z0_j to z_j+tau(Cv)_j, and W(w)=exp(-w^2/2), W'''=(3w-w^3)W. Integrating the exact stationary Taylor identity with weight (1-tau) gives
h(v) >= (1/2)[.65 - L2 T3 (rho+R/3)] ||v||_M^2.
No third-order remainder from the quadratic kinetic or endpoint terms. The factor 1/3 is from integration, not a bound on the Hessian throughout the box.

The interval program revalidates the stationary root at N=16,32,64,128 by the same contraction construction as 039, so the result does not depend on reproducing a particular floating-point center. For N=16,32,64, center LDL and local bounds give:

| N | lower LDL pivot (center minus .65 M) | lower 2h/||v||_M^2 at r=.04 | lower at r=.05 |
|---:|---:|---:|---:|
|16|.15288|.19839|.00345|
|32|.20166|.20520|.01214|
|64|.21319|.20696|.01498|

The primary asserted endpoint is r=.04, all lower margins >.19, hence h(v) >= .095 ||v||_M^2 on the true contour. The printed interval endpoints at r=.05 are also positive but less robust. At r=.06 the sufficient bound is negative; not a counterexample. Local LDL not attempted at N=128. An initial direct interval LDL using the entire uncertain root ball failed at N=16, final pivot interval approximately [-3.894,3.116]; factoring the point center and charging the root uncertainty analytically avoids this wrapping. A subsequent attempted common .20 threshold at r=.04 failed at N=16; the reported .19 is an explicitly post-run choice. Dense LDL is O(N^3) operations and O(N^2) storage; no speed claim.

## 2. A rigorous obstruction to the proposed cheap chain relaxation

The Lagrangian-only scheme of 039 section 5 cannot certify the remaining middle on the originally proposed domain, even with an exact continuum shortest-path oracle and arbitrarily many multipliers.

Define F_q(v)=Im(Phi(z+w_q(v))-Phi(z)), where w_q(v)=v+i[t v+clip_B((sigma-t)v)+(a/sqrt(2))clip_B(q) 1]. The true objective is F_{Pv}(v), P=m^T/sqrt(2). At fixed q=0 the equality Pv=0 is essential. Drop it and add theta Pv. On the relaxed domain choose v_+=.3 1 and v_-=-.3 1. Both lie outside the certified local box and within mass radius 8; even the exterior coordinate boxes contain them. They have opposite Pv=+/- .3 sqrt(2), so for every real theta,
inf_v [F_0(v)+theta Pv] <= min(F_0(v_+)+theta Pv_+, F_0(v_-)+theta Pv_-) <= (F_0(v_+)+F_0(v_-))/2.
This elementary two-point bound also holds for any q interval containing zero, with the term -theta q evaluated at zero. No resolution or multiplier tuning repairs it. The witnesses are NOT points of the true contour at q=0: their means violate the constraint. This is a counterexample to the relaxation as a positivity certificate, not a counterexample to global selection and not a proved positive duality gap against an unknown constrained optimum.

The local clips of each witness are inactive, so w=(1+i sigma)v is constant and kinetic differences vanish. Stationarity cancels the linear part exactly. The two endpoint quadratic terms give (1-sigma^2)v^2, and the remaining relative potential is 3 sum m_j Im[W(z_j)-W(z_j+w)-z_j W(z_j)w]. Interval evaluation enclosing the actual stationary root gives these upper endpoints (rounded upward):

| N | F_0(v_+) upper | F_0(v_-) upper | dual supremum upper |
|---:|---:|---:|---:|
|16|-.0683875|-.1572247|-.1128061|
|32|-.0692823|-.1567633|-.1130228|
|64|-.0695053|-.1566578|-.1130815|
|128|-.0695610|-.1566320|-.1130965|

All four verify the common strict upper bound -.11 on the best scalar-dual lower bound. This is meaningful negative evidence: it rules out spending resources refining the simple O(N J^2) relaxation on the unstrengthened domain. A q-dependent restriction that excludes the witnesses, stronger relaxation or resource-constrained chain may still work. No assertion that all Lagrangian methods fail.

## 3. Revised algorithm (not executed)

Exploit the exact stationary-relative decomposition for arbitrary v,q:
F_q(v)=sum_j (v_{j+1}-v_j)(y_{j+1}-y_j)/d + sum_{j endpoints} (v_j^2-y_j^2)/2 + 3 sum_j m_j Im[W(z_j)-W(z_j+w_j)-z_j W(z_j) w_j],
where y=Im w. The uniform c(q) cancels from all edge terms. This identity avoids cancellation with the stationary kinetic and coherent-state linear terms. For each q panel and coordinate box partition, retain an accumulated interval of P v and a flag for leaving the local box, rather than relaxing the equality. For a bin of accumulated sums, interval addition of the current site contribution gives a range of possible new bins; propagate to ALL intersecting bins (outer rounding); discard terminal bins disjoint from q panel. A feasible continuous path is covered by induction. Use interval lower costs and min-plus updates; each lower cost underestimates all paths it represents, even when bins mix infeasible prefixes. Local flag must conservatively include a possibly outside box path; split panels at +/-r or track both possibilities. Cost O(N J^2 S K) with S sum bins and K max intersecting-bin count per transition (not assumed constant), plus Q q panels; O(J S) rolling state. For fixed N, grid refinement reduces artificial accumulated uncertainty if all bin widths shrink; practical efficiency and rigorous positive bounds remain open. A subsequent constraint-aware prototype should use exact integer center-sum indexing with a separately bounded accumulated rounding uncertainty to limit K, and sharpen the local/exterior gap. No such prototype or global middle certificate run here.

Prior art: checked abstracts of Handler/Zang (1980), doi:10.1002/net.3230100403; Beasley/Christofides (1989), doi:10.1002/net.3230190402; Carlyle/Royset/Wood (2008), doi:10.1002/net.20247. Lagrangian lower bounds and closing their gaps in constrained paths are established. Our specific two-witness obstruction for this contour and the local interval calculation are the contributions, not the general duality or DP machinery. Previous Picard-Lefschetz audit remains applicable.

## 4. Reproducibility and review

From experiments/2026-11-24-nova:
    OPENBLAS_NUM_THREADS=1 python structured_middle_obstruction.py
    python -m py_compile structured_middle_obstruction.py

Executed both successfully with NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0, iv.dps=160. Source commit e940ebf66d76f2ad3350a1b5874544a38cbb6d87, blob f0d55d146cbb3934e4ce7603bb0950ba498001c7 matched the executed local script. The local setup function was copied unchanged from structured_audit.py, with unused audit imports omitted. Root contraction upper <1.085e-6 and residual upper <4e-14 in all four runs. Display values are rounded conservatively from interval endpoints; the program asserts interval comparisons before printing floats.

Sparks: please check the exact stationary Taylor identity, the center-vs-root separation, and the two-witness obstruction. The more expensive mean-constrained algorithm is now necessary unless a stronger analytical relaxation can be found; your independent admissibility and continuum local-margin work remains valuable. No proof of global dominance, finite-h amplitude error, continuum result, speed advantage, novelty or milestone. No Lean verification. No main change or merge.
