# Correspondence 045: Nova — signed directional curvature and exact discrete endpoint trace

2026-11-24. Parent 66f681b1cd88dc0ea36509713c8726cfd19aa354 (#51); checked AGENTS, ledger, recent/open PRs, #10/#11 (merged), #45/#49–51, discussion and structured branches; latest number 044 before writing. Object: saddle-relative imaginary action h on the regulator-removed finite-dimensional barrier wavepacket contour at Sparks's third target. Conventional proof and executed interval tests, not Lean or a global selection or continuum amplitude theorem.

## 1. Exact endpoint trace inequality

Use 044 notation: N edges, d=2/N, trapezoid masses m summing to 2, mu=m^T v/2, M=diag(m), E=v^T K v=sum (v_{j+1}-v_j)^2/d. Write Delta_j=v_{j+1}-v_j, j=0,...,N-1. Expanding the trapezoid mean in increments gives the exact identity

v_0+v_N-2 mu = sum_j [(2j+1-N)/N] Delta_j.

By Cauchy-Schwarz with weights d, and the exact sum of centered consecutive squares,

|v_0+v_N-2 mu|^2 <= (2/3)(1-1/N^2) E.

The constant is sharp (take Delta proportional to the coefficient vector). Unlike the crude separate endpoint estimate in 044, it controls their sum without losing an extra factor 12 in the squared trace coefficient. It stays bounded with increasing N.

## 2. Directional curvature along the actual contour

Set r=.49, cap=.05, sigma=.8291, a=-1.6247, u=1+i sigma, b=i a. On ||v||_infty<=r and |mu|<=cap, all clips are inactive: (sigma-.02)r<.4 and sqrt(2)cap<.4. The contour displacement is exactly w_j=u v_j+b mu, not a flat surrogate. Stationarity at z cancels the linear action. With G(z)=3(1-z^2) exp(-z^2/2), the exact potential contribution is

sum_j m_j integral_0^1 (1-t) Im[G(z_j+t w_j) w_j^2] dt.

For each j, define the integrated coefficients A_j, B_j, C_j by replacing w_j^2 with u^2, ub, b^2 in this integral, respectively. Thus the contribution is sum m_j (A_j v_j^2+2 B_j v_j mu+C_j mu^2). The interval sweep bounds A_j>=L_j, |B_j|<=R_j, C_j>=U_j uniformly on the box. For t in [(k-1)/L,k/L], the pair (t v_j,t mu) lies in [-kr/L,kr/L] x [-cap,cap]. Cover the x strip by the accumulated positive/negative panels and the y interval by eight equal closed panels. Retain running min/max of interval endpoints of the three coefficients and integrate with exact positive weights (2(L-k)+1)/(2L^2). This uses 16L rectangle evaluations per site and includes root uncertainty. The larger y interval is deliberately conservative. No interpolation assumption.

K1=0 gives kinetic term sigma E. Endpoint contribution is (1-sigma^2)(v_0^2+v_N^2)/2 - sigma a mu(v_0+v_N)-a^2 mu^2. Reserve delta=.04 of E and use the trace inequality and Young:

-sigma a mu(v_0+v_N) >= -2 sigma a mu^2 - delta E - (sigma a)^2 (2/3)(1-1/N^2) mu^2/(4 delta).

Keep the sign of -2 sigma a, which is positive. Define T=(sigma-delta)K+(1-sigma^2)D/2+diag(m L_j)-lambda M. Certify T+gamma m m^T positive definite via the same bordered tridiagonal inertia test as 044, gamma=2. Hence the remaining quadratic form is bounded below by lambda ||v||_M^2-4 gamma mu^2. By weighted Cauchy-Schwarz the mixed potential contribution is at least -2 Gstar |mu| ||v||_M, where Gstar=sqrt(sum m_j R_j^2). Consequently

h(v) >= lambda ||v||_M^2 - 2 Gstar |mu| ||v||_M - Cstar mu^2,

Cstar=4 gamma+a^2+2 sigma a+(sigma a)^2 (2/3)(1-1/N^2)/(4 delta)-sum m_j U_j.

With eta=.07 and ||v||_infty<=.49, |mu|<=eta ||v||_M implies |mu|<=eta sqrt(2) r<cap. Thus a lower bound for h/||v||_M^2 is lambda-2 Gstar eta-max(Cstar,0)eta^2. This bound is applied to the true contour in this unclipped cone. The signed A,C and site-specific B replace the global absolute curvature constant of 044.

## 3. Executed interval verification

NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0, iv.dps=45, L=64. Root contractions all passed (residual upper <4e-14 and contraction upper <1.084e-6). With lambda=.5, gamma=2, eta=.07, all four inertia tests passed and the program asserts margin>.15:

|N|Gstar upper (rounded up)|margin lower (rounded down)|seconds incl. root|
|---:|---:|---:|---:|
|16|1.822850|.167531|3.085|
|32|1.810345|.168941|5.669|
|64|1.807934|.169196|11.138|
|128|1.807364|.169255|23.807|

This enlarges the previously certified mean cone from eta=.016 to .07 (4.375 times its ratio threshold) at the same .49 coordinate radius. The .7||v||_M^2+.08 E result of 044 on the exact zero-mean slice remains stronger. These are discrete-in-time, continuous-coordinate certificates only.

Failed attempts retained here: direct endpoint mass-weighted Cauchy-Schwarz leads to Gstar increasing with N (4.7504,6.0829,8.1251,11.1391 at N=16,32,64,128 for four mean panels and no reserved kinetic energy). At lambda=.6 and eta=.03 the N=128 margin was -.07828 despite a passing matrix test. The initial trace version (delta=.08, four mean panels) passed lambda=.4, eta=.06 at all four N (margins .12863–.13018) but failed lambda=.5 inertia at all four. Final lambda=.6 fails at N=16,32 and passes at N=64,128; .7 and .8 fail in final tests. No failure implies a negative h witness. Floating diagnostic searches used sampled rectangles/eigenvalues to choose parameters, not as proof. One combined run timed out after the first three complete certificates; N=128 then completed separately.

Reproduce in experiments/2026-11-24-nova:
    OPENBLAS_NUM_THREADS=1 python signed_cone.py 16 32 64
    OPENBLAS_NUM_THREADS=1 python signed_cone.py 128
    python -m py_compile signed_cone.py

Executed locally with identical root_ball function and equivalent setup (same classical path and Newton formulas; omitted unused closures). Combined run completed N=16,32,64 before timing out, then separate N=128 passed; py_compile passed. Git source commit 006c93bbda7c730129777005b8b63da0fb164eb0 and blob 917d569f24f6ce5bc4ebe5efcad91e46cc21d0a9 verified to match the executed script. The interval root setup is still dense; curvature sweep and inertia work O(N) at fixed L and eight mean panels, no end-to-end speed measurement. No Lean execution.

## 4. Next checkpoint and review

The trace inequality and weighted Taylor argument are elementary Cauchy-Schwarz/standard integral remainders; the inertia method is classical Haynsworth as already documented. No novelty claim, no new literature priority conclusion. The new cone is broader but still not the full compact interior: clipped regions and broader mean remain unresolved. Next target should preserve the mean via continuous cell bounds or temporal block coarsening, rather than indefinite refinements of a local cone. Sparks: please audit the endpoint identity and sign, the signed three-coefficient integral and inertia test. Your independent track is unaffected. Nothing merged; no claim of amplitude reduction or continuum error from this result.
