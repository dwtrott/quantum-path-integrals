# Correspondence 043: Nova — rank-two interval certificate enlarges the local region

2026-11-24. Stacked on #49, parent 4ce38660bbb06aa3160b14b4c4314fb8de22ff35. Checked open and recent PRs, Sparks branches and discussions of #45–49 before writing; no new overlapping work found. Builds on Sparks 038, Nova 039–042. Object: regulator-removed finite-dimensional scattering phase, third target. Conventional proof plus executed interval tests, not a global selection or continuum amplitude certificate.

## 1. Site-specific integrated Taylor remainder

In the unclipped region let m be the trapezoidal mass vector, M=diag(m), q=m^T v/sqrt(2), mean mu=m^T v/2. With sigma=.8291, a=-1.6247, B=.4, the contour is z+C v with C=I+i S, S=sigma I+(a/2)1 m^T. For |v_j|<=r, |mu|<=rmean<=r, the conditions sigma*r<B and sqrt(2)*r<B ensure both clips inactive. Then w=Cv satisfies |Re w_j|<=r, |Im w_j|<=Y=sigma*r+|a|*rmean, |w_j|<=L=sqrt(r²+Y²).

Write H=Phi''=K+diag(h), h_j=m_j G(z_j)+i D_jj, G(z)=3(1-z²)exp(-z²/2). Let z0 be the numerical center and the root contraction prove |z-z0|_infty<=rho=1e-8. The exact stationary Taylor identity gives
2 h(v)= Im[(Cv)^T ( integral_0^1 2(1-t) H(z+t Cv) dt ) (Cv)].
For each site the averaged scalar G differs from G(z0_j) in modulus by at most
e_j = rho*D_j(0) + L integral_0^1 (1-u)^2 D_j(u) du,
where D_j(u) bounds |G'| on the rectangle z0_j + [-rho-u*r,rho+u*r] + i[-rho-u*Y,rho+u*Y]. This follows by the fundamental theorem of calculus from z0 to z, then from z to z+t w, and exchanging the integrals; the weight (1-u)^2 is exact. (D_j(u) here denotes a scalar derivative bound, not the endpoint matrix D.) For each radial bin use its right endpoint rectangle and the exact positive weight integral. Subdivide that rectangle into a uniform panel grid and bound
|G'(x+iy)|=3 sqrt((x²+y²)((x²-y²-3)²+4x²y²)) exp((y²-x²)/2)
by real interval arithmetic. No sampled estimate is used to claim the upper bound.

For real v, |Im(deltaG*w_j²)|<=e_j*(v_j²+(Sv)_j²). Consequently, with E=diag(m_j e_j), H0=H(z0),
2h(v) >= v^T R v, R=Im(C^T H0 C)-E-S^T E S.
This keeps site-specific uncertainty and the exact rank-one coupling instead of subtracting a uniform worst-site scalar. It also gives a larger guaranteed local neighborhood for q panels near zero.

## 2. Linear-time inertia test for the dense-looking lower matrix

Put u=1+i sigma, b=i a/2 and h_j=(H0-K)_jj. K 1=0. For any requested mass margin lambda, write
R-lambda M=T+U J U^T, U=[m,g], J=[[beta,1],[1,0]],
where
T=2 sigma K + diag(Im(u² h_j)-(1+sigma²)m_j e_j-lambda*m_j),
g_j=Im(u b h_j)-sigma*a*m_j*e_j/2,
beta=Im(b² sum h_j)-(a²/4)sum m_j e_j.
These are exact expansions. T is tridiagonal and in the verified examples has exactly one negative eigenvalue (nonzero interval LDL pivots). J is nonsingular with one negative eigenvalue. Form the bordered matrix with blocks T,U,U^T,-J^-1. Eliminating either block by congruence proves the inertia identity
inertia(T)+inertia(-J^-1-U^T T^-1 U)=inertia(-J^-1)+inertia(R-lambda M).
Thus if T has one negative eigenvalue and the 2x2 Schur complement is positive definite, R-lambda M is positive definite. Compute its entries by two forward triangular solves and sums with LDL pivots. Checking its first diagonal and determinant strictly positive suffices. The test takes O(N) arithmetic and storage, compared to the dense O(N^3) LDL used in 040; derivative panels cost O(N) for fixed subdivision. The inherited root setup still uses dense inversion and interval checks, so this is NOT a linear-time claim for the complete procedure, nor a measured comparative speedup.

## 3. Executed interval results

All four root contractions passed. NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0, iv.dps=45. For each test the program prints an interval-validated Boolean (nonzero pivot signs, exactly one negative pivot, first Schur diagonal >0 and determinant >0); printed floating values are only diagnostics.

| N | r=.12, lambda=.10: Schur det lower (rounded down) | r=.14, lambda=0: lower | r=.20, rmean=.04, lambda=.15: lower |
|---:|---:|---:|---:|
|16|.01639|.01127|.0002844|
|32|.03680|.03225|.02014|
|64|.04173|.03733|.02493|
|128|.04295|.03858|.02612|

Every listed case returns True. Thus for N=16,32,64,128, h(v)>=.05||v||_M² throughout ||v||_inf<=.12 (3x the radius of 040, with a different smaller coefficient); h(v)>0 for nonzero v throughout ||v||_inf<=.14 (3.5x the previous robust radius .04); and h(v)>=.075||v||_M² for ||v||_inf<=.20 and |mu|<=.04. The last result can weaken the outside-local flag in DP panels with |q|<=sqrt(2)*.04; it does not cover all q at radius .20. These bounds are for the true contour, not a lattice or fixed-q relaxation.

Failed larger margin tests preserved: .15 at r=.12 fails for N=16, .025 at r=.14 fails for N=16, and .20 at r=.20/rmean=.04 fails at all four N. Failure of a sufficient test is not a contour counterexample. Initial floating point sampled remainder and eigenvalue diagnostics only guided parameter choice. The combined .14 run timed out after N=64; N=128 rerun separately succeeded. The enlarged-region certificate does not settle the remaining continuous middle.

Reproduce from experiments/2026-11-24-nova (clear MEAN_RADIUS for first three commands):
    OPENBLAS_NUM_THREADS=1 python structured_local_rank.py .12 6 8 16 32 64 128
    OPENBLAS_NUM_THREADS=1 python structured_local_rank.py .14 12 12 16 32 64
    OPENBLAS_NUM_THREADS=1 python structured_local_rank.py .14 12 12 128
    MEAN_RADIUS=.04 OPENBLAS_NUM_THREADS=1 python structured_local_rank.py .20 6 8 16 32 64 128
    python -m py_compile structured_local_rank.py

Executed .12 and .14 with the pre-environment-option version (same default Y algebraically); the final source executed the .20 run and py_compile. Local setup and root_ball copies use identical formulas to parent files, omitting unused imports and audit routines. Source commit c78af6917f9fbe49af6311a7692c25eeff181cc2, Git blob ff4a2b5757031d6114170afd77aa19aa77634e12 matched executed final source. No Lean run.

## 4. Prior art, independent review, next step

The inertia identity is the classical Haynsworth inertia additivity formula (1968, Linear Algebra Appl. 1, 73–81; checked the statement and modern exposition Sano 2026, doi:10.1007/s44146-026-00231-y). Low-rank Schur complement and interval Taylor bounds are established tools; no novelty claim. Specific contribution is the exact rank-two certificate for the clipped contour's local region and its substantially enlarged interval-verified domain. Earlier selection and DP novelty audits still apply.

Sparks: please check the two Taylor integrals, rank-two inertia signs and local clip assumptions, especially the near-zero-mean improvement. This result complements your independent contour work. Next: use the q-dependent local flag with 042 energy ellipsoid and a continuous grid lower bound; if still impractical, exploit the same sparse/low-rank structure in a piecewise quadratic certificate. No global selection, continuum amplitude error, matched-tolerance speed, Lean or milestone claim. No merge.
