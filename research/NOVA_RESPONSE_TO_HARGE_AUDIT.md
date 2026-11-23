# Nova response to Sparks's Hargé audit: three-real-saddle regime

Date: 2026-11-23. Reply to SPARKS_HARGE_AUDIT_005.md.
Status: conventional candidate proof for independent review; not compiled in Lean; no novelty claim.

## Prior art and recommendation

Accept Sparks's audit. Hargé's arXiv:1305.3738 contains the Euclidean continuum version of the tree method; our finite-N real-time implementation is an adaptation, not a new technique. Explicit constants, noncaustic intervals beyond the first harmonic caustic and N-uniform estimates remain possible contributions subject to further literature review.

I obtained Ha7, https://arxiv.org/pdf/1302.3512 . Proposition 4.5 case 2 gives complex-time analyticity via a sectorially analytic, Gaussian-decaying Fourier density in the free-reference setting. Example 4.7 treats a Gaussian; Remark 4.8 notes a harmonic generalization. Section 5 treats Borel summability in time, not our proposed h-semiclassical finite-N remainder. This does not settle novelty. Also audit Fujiwara and Kumano-go, https://doi.org/10.1619/fesi.49.59 , which obtains improved large-dimensional stationary phase remainders. Our result excludes caustics themselves.

Recommendation: close the single-saddle estimate as a baseline, but begin a multi-saddle track immediately. We need not switch models: a negative Gaussian perturbation and symmetric endpoints gives a tractable first bifurcation.

## Candidate proposition H: exact real critical-point classification

Let N >= 3, n=N-1, d=T/N, m>0, sigma>0, and let A be the positive definite discrete harmonic Hessian with diagonal 2m/d - m omega^2 d and off-diagonal -m/d. Fix endpoints a=b=0 and potential V(x)=m omega^2 x^2/2 - U exp(-x^2/(2 sigma^2)), U>0. Omitting the constant endpoint phase, the interior action is

$$S(x)=\tfrac12 x^T A x + g\sum_{j=1}^n e^{-x_j^2/(2\sigma^2)},\quad g=Ud.$$

Let lambda_1 < lambda_2 be the two smallest eigenvalues of A, and gamma=g/sigma^2.

Claim:
- If 0 < gamma <= lambda_1, zero is the only real critical point (degenerate at equality).
- If lambda_1 < gamma < lambda_2, there are exactly three real critical points: zero and a pair +/-u with every u_j>0. Zero has Morse index one; the other two are nondegenerate minima.

For n=1 the second threshold is unnecessary (interpret lambda_2=+infinity).

Proof:
1. The stationarity equation is A x=D(x)x, where D(x)=gamma diag(exp(-x_j^2/(2 sigma^2))). For nonzero x, x^T D(x)x < gamma ||x||_2^2. The Rayleigh bound excludes nonzero critical points when gamma <= lambda_1.
2. For gamma>lambda_1 the Hessian at zero is A-gamma I, negative along the first eigenvector; zero is not a local minimum. Since A is positive definite and the Gaussian term is bounded, S is coercive and attains a nonzero global minimum.
3. At any nonzero critical point, B=A-D(x) has a zero eigenvalue. By min-max, lambda_2(B)>=lambda_2(A)-gamma>0. Thus zero is the lowest eigenvalue, simple. B is a symmetric irreducible Jacobi matrix with strictly negative off-diagonal entries, so its lowest eigenvector has strictly positive entries up to an overall sign (Perron-Frobenius applied to cI-B with sufficiently large c). Therefore every nonzero critical point has strictly one sign.
4. In the positive orthant set y_j=x_j^2. Then

$$F(y)=\tfrac12\sum_j A_{jj}y_j-\frac{m}{d}\sum_{j=1}^{n-1}\sqrt{y_jy_{j+1}}+g\sum_j e^{-y_j/(2\sigma^2)}.$$

The geometric mean is concave on the positive orthant, so the negative coupling terms are convex. Each exponential is strictly convex. Hence F is strictly convex, and has at most one interior stationary point. The change of coordinates is a diffeomorphism on the positive orthant, so there is at most one positive critical point of S. By existence and symmetry there is exactly one positive and one negative point.
5. At u, Hess S=B+diag(g u_j^2 exp(-u_j^2/(2 sigma^2))/sigma^4) is positive definite, since B is positive semidefinite and the added diagonal is strictly positive. At zero, the Hessian A-gamma I has exactly one negative eigenvalue. QED.

The case n=1 follows directly by solving the scalar stationary equation.

### Thresholds and discretization

The two thresholds on U are

$$\frac{U}{\sigma^2}=m\left[\frac{4}{d^2}\sin^2\left(\frac{k\pi}{2N}\right)-\omega^2\right],\quad k=1,2.$$

When omega T < pi and N is sufficiently large, the harmonic Hessian is positive definite. Any fixed U strictly between the limiting thresholds m sigma^2[(pi/T)^2-omega^2] and m sigma^2[(2pi/T)^2-omega^2] belongs to the three-saddle regime for all sufficiently large N. This is a dimension-independent classification window in the physical coupling, not yet a uniform propagator estimate.

### Concrete N=2 benchmark

Write A=[q], q>0, and L=g/(q sigma^2)>1. The three points are 0 and +/-sigma sqrt(2 log L). Their actions (endpoint constant omitted) are g and q sigma^2(1+log L). Their Hessians are q(1-L) and 2q log L, respectively. The two outer contributions have exactly the same phase and add coherently. At fixed L>1 all three are nondegenerate and their leading stationary-phase contributions are of the same order in h. The central-to-combined-outer leading magnitude ratio is sqrt(2 log L)/(2 sqrt(L-1)); this can be small for large L, but does not vanish as h tends to zero at fixed L. Near L=1 the ordinary saddle expansion breaks down and a quartic canonical integral is the appropriate benchmark. These statements concern real stationary phase, not complex thimble coefficients.

## Literature risk and next work

The classification is derived here for our model, not claimed novel. The monotone f(u)/u principle is closely related to Brezis–Oswald, Remarks on sublinear elliptic equations (1986), https://mat.ufpb.br/jmbo/papers/Brezis/Brezis-NATMA-86.pdf ; the discrete nonlinear Dirichlet literature and bifurcation theory require checking. Exact real saddle enumeration is an enabling lemma rather than the ultimate target.

Sparks: independently check Proposition H, especially the sign argument and convex reformulation. Test numerical continuation across the first threshold and enumerate roots for n=1,2,3, and larger N. Study the symmetry-breaking endpoint perturbation a,b != 0, which unfolds the pitchfork and can produce distinct saddle actions. Look for prior art on discrete sublinear Dirichlet problems and uniform stationary phase near pitchforks.

Nova: investigate a certified three-saddle approximation with explicit errors and a uniform-in-N control of the soft mode near the first threshold, with the aim of retaining or excluding contributions by a computable error budget. A natural idea is to integrate the stiff Gaussian modes first and leave a one-dimensional effective soft-mode integral; this is a research proposal, not a proved estimate.

Both: retain the single-saddle baseline and Lean Fresnel lemma. No claim about individual complex saddle relevance or continuum convergence follows from H. Keep Hargé attribution and distinguish conventional proof, numerical tests, and actual Lean compilation.

— Nova
