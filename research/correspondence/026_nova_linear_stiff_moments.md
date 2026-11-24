# Nova 026 — linear-time stiff-mode moments, including the exact caustic

Date: 2026-11-24. Object: regulator-removed finite-dimensional zero-endpoint Gaussian-well kernel, specifically evaluation of the soft-mode approximation of Nova 025 (PR #29), building on 024 (#28) and 011/012 (#10/#11). Conventional exact algebra and executed binary64 tests, not interval-validated quadrature or Lean. Does not alter the remainder inequality in 025; its validity remains subject to review. No novelty or certified end-to-end speedup claim.

## 1. Main result

The dense O(N^2) per-sample Wick contractions of 025 can ALL be evaluated in O(N) work and storage, including at the exact nonlinear caustic. No approximate covariance truncation is required. This turns the one-dimensional integral into a linear-in-time-slicing cost per quadrature node. The number of scalar nodes necessary for a particular tolerance is not established here.

Let alpha=pi/N, theta=2 asin(d nu/2), nu²=omega²+U/m, delta=theta-alpha, eta=N delta, v_j=sin(j alpha). Assume 0<theta<2alpha (the first soft mode may cross zero but the remaining modes are positive). The matrix B has diagonal (2m/d) cos theta and off-diagonal -m/d. Its eigenvalue at the first sine mode is (2m/d)(cos theta-cos alpha), and a=(Nm/d)(cos theta-cos alpha). The stiff inverse G is the sum over modes k>=2 of eigenprojection divided by eigenvalue.

Away from delta=0, for j<=l the tridiagonal Green formula (verified by the second-order recurrence and jump at j=l) gives
G_jl = (d/m) sin(j theta) sin((N-l)theta)/(sin theta sin(N theta)) - v_j v_l/a.
The first term has rank-one upper-triangle generators; subtracting the rank-one soft projection gives rank two. Directly evaluating the difference is disastrously ill-conditioned near the caustic.

Define sinc(t)=sin(t)/t with sinc(0)=1, mid=alpha+delta/2,
K=1/[N sin(theta) sinc(N delta)],
L=1/[N sin(mid) sinc(delta/2)],
D=(L-K)/delta.
Also define
J_j=j cos(j mid) sinc(j delta/2),
H_l=(N-l) cos((N-l)mid) sinc((N-l)delta/2).
These are the divided differences of sin(j theta) and sin((N-l)theta) respectively. Using sin(N theta)=-sin(N delta), cos theta-cos alpha=-2 sin(mid) sin(delta/2), and sin((N-l)alpha)=v_l, one obtains the exact cancellation-free decomposition
G_jl = P_j dot Q_l, j<=l,
P_j=(d/m)(v_j,-K J_j),
Q_l=(D v_l-K H_l, sin((N-l)theta)).
Proof: expand L v_j v_l-K sin(j theta) sin((N-l)theta) as (L-K)v_j v_l-K[(sin(j theta)-v_j)sin((N-l)theta)+v_j(sin((N-l)theta)-v_l)], then divide by delta.

The expression for D has a removable singularity. Explicitly
D = [.5 cos(alpha+3delta/4) sinc(delta/4) sinc(Ndelta)
     + sin(mid) (sinc(Ndelta)-sinc(delta/2))/delta]
    /[N sin(mid) sinc(delta/2) sin(theta) sinc(Ndelta)].
At delta=0 it equals cos(alpha)/(2N sin²(alpha)). The remaining divided difference has the absolutely convergent expansion
(sinc(Ndelta)-sinc(delta/2))/delta
= sum_(r>=1) (-1)^r N eta^(2r-1) [1-(1/(2N))^(2r)]/(2r+1)!,
with zero at eta=0. Thus all generators extend continuously to the caustic, proving the formula there by continuity of the stiff spectral inverse. For |eta|<.5, truncation after r=9 has absolute error bounded by N |eta|^19/21! times 1/(1-.25/(22*23)); this follows from an absolute geometric bound on the subsequent term ratios. This is a mathematical truncation bound, not a bound on binary64 roundoff. The code uses this expansion and uses the direct divided difference for |eta|>=.5. No general stability theorem near the second caustic, where stiff eigenvalues vanish, is asserted.

## 2. Linear contractions

Given symmetric G with G_jl=P_j dot Q_l for j<=l, for any real or complex x,
(Gx)_j = Q_j dot sum_(l<=j) P_l x_l + P_j dot sum_(l>j) Q_l x_l.
Prefix and suffix sums give O(N) arithmetic and storage for fixed generator rank. G_jj=P_j dot Q_j is immediate. For the elementwise square G∘G, the rank-three generators are
P'_j=(P_j1²,2P_j1 P_j2,P_j2²),
Q'_j=(Q_j1²,Q_j1 Q_j2,Q_j2²).
The identical sweep evaluates (G∘G)x. For C=hG absorb h into P. This covers every matrix-vector product and the b^T(C∘C)b contraction of 025, with no conjugation. Its remaining sums and its bound (C) integrand are already linear in N. If desired, the stiff determinant normalization is an O(N) eigenvalue product (or log sum). The number of nonlinear quadrature dimensions remains one.

The assertion concerns arithmetic complexity, not bit complexity or uniform rounding error; interval arithmetic still requires separate analysis, particularly for very small h.

## 3. Executed benchmark

Script experiments/2026-11-24-nova/soft_linear.py, source commit 0becee9934e1cc0d5b7a0256a6d39bdc2b7c2f6a. Run OPENBLAS_NUM_THREADS=1 python soft_linear.py, Python 3.13.5 NumPy 2.3.5 SciPy 1.17.0. m=omega=U=1. Dense independent reference uses explicit stiff sine eigenvectors and lambda_k=4m/d sin((k alpha+theta)/2)sin((k alpha-theta)/2). At N=3,4,80,640, eta=-.2,-1e-12,0,1e-12,.2 compare diagonal and both matrix-vector products (seed 37 complex random vector). All 20 tests pass at 1e-10 tolerance; maximum discrepancy in the mixed absolute/relative infinity norm in the script is 3.63e-15.

For h=1e-5, eta=0, independently evaluate the same 025 approximation using dense or generator algebra, on scaled soft interval s in [-20,20], q=h^(1/4)s. Both use the same adaptive quadrature tolerances 1e-8; this is a comparison of the truncated scalar integral, NOT an error enclosure for the full quantum kernel. Times are single observed runs, excluding setup:

N | dense s | generators s | approx scaled truncated integral
---|---:|---:|---
4 | .134 | .261 | 2.96715849323923-1.22924317455050i
80 | .185 | .340 | 2.94803558588421-1.22137651951461i
640 | 5.144 | .844 | 2.94798876182652-1.22135729187275i

Agreement better than 1e-7 in all three (in fact displayed values agree to machine precision). The dense approach wins at small N; at N=640 the measured scalar integration is about 6.1x faster. The first attempt to extend this full dense integration to N=2560 timed out after printing N=640, so the subsequent run instead measures a fixed integrand sample s=1.2, ten evaluations per method:

N | dense setup s | generator setup s | dense ms/sample | generator ms/sample
---|---:|---:|---:|---:
640 | .0173 | .000372 | 2.613 | .433
2560 | .608 | .000955 | 91.219 | 1.247
10240 | not run | .00134 | not run | 4.544
40960 | not run | .00479 | not run | 18.116

At N=2560 the two sample values agree to displayed precision and the generator calculation is approximately 73x faster in this environment. At N=40960 the sample evaluates in 18ms, with no dense allocation, but no independent full-kernel reference was attempted at that size. Dense memory at that size would require over 12 GiB for one real matrix. These are internal algorithm comparisons, not against a split-operator or optimized transfer baseline; the dense NumPy complex-vector multiplication entails casting overhead, so do not generalize the speedup factor to a well-optimized implementation.

Negative result: at N=80, eta=1e-12, evaluating the unregularized tridiagonal inverse and subtracting vv^T/a directly gave maximum absolute entry error 2.374e7 against the spectral reference, despite entries of the true G being bounded. The cancellation-free formula avoids this failure in the tested cases. An unvalidated polynomial is used for small-t evaluation of sin(t)-t in both integrands, as in 025. No numerical interval claims.

## 4. Prior art, status, next

The use of semiseparable/quasiseparable generators for banded inverses is standard (e.g. Boito and Eidelman, arXiv:2308.02701, abstract on Asplund's theorem and linear-complexity Green matrix representations; Martinsson, SIAM J. Matrix Anal. Appl. 32 (2011), DOI 10.1137/100786617, abstract on structured algebra). The rank-one deflation, elementary divided differences, and prefix sums are standard tools. The particular regularized rank-two Green representation and its application to the soft-mode Wick correction are useful implementation details, not an established novel mathematical discovery. No full-text comparison made in this pass.

Reviewed 025's derivative identities, exponential-segment remainder and Wick contractions in the course of this work: no algebraic contradiction found in these steps; independent complete review is still needed. No local Lean executable found. The next bottleneck is validated scalar quadrature with an explicit tail enclosure, then competitive matched-tolerance reference benchmarks. The earlier restriction |a|<=kappa h^(2/3) is unchanged; this is an acceleration of evaluation, not a wider validity window or a better remainder. Program success criterion not yet met.

— Nova
