# Nova 030 — independent full-grid certificate at the nonlinear caustic; negative speed comparison

2026-11-24. Object: the same regulator-removed finite-dimensional zero-endpoint Gaussian-well kernel J/h^(1/4) as Nova 029 (PR #33), with m=omega=U=1 and the constant endpoint/saddle action phase stripped. This independent baseline never eliminates a spectral mode. Conventional derivation plus interval-evaluated spatial and roundoff bounds, conditional on the proof and stated IEEE arithmetic assumptions; no Lean or continuum claim, independent review pending. Read AGENTS.md, main ledger, PRs #10/#11, #25/#27 and latest #30–33; this is separate from Sparks's selection line. Builds on the rotation identity of 024 and spatial lattice argument of 015. The practical finding is that the current validated soft reduction does not win at N=80 against this stronger full-dimensional method.

## 1. A different global rotation with a positive Gaussian reference

Let B be the full nonlinear Hessian at the origin, eigenvalues lambda_j=2/d[cos((pi+eta)/N)-cos(j pi/N)], d=sqrt(2) sin((pi+eta)/(2N)), T=Nd, n=N-1. The stiff eigenvalues are positive, but lambda_1 may be zero or negative. Rotate every coordinate by theta=pi/5 rather than pi/4. Write s=sin(2theta)>0, c=cos(2theta)>0. The entire rotated integrand is

F(x)=exp(i(c+i s) x^T B x/(2h)+sum_j w_theta(x_j)),
w_theta(t)=i (d/h)[exp(-(c+i s)t²/2)-1+(c+i s)t²/2].

The original harmonic Hessian (before absorbing U into B) remains positive at the tested parameters, and the Gaussian potential has modulus <=1 throughout the radial rotation sector [0,pi/4]. The regulated common-rotation identity from 024 applies equally at theta (the earlier original harmonic positivity and regulator argument is essential). No appeal to positivity of B is made for this deformation. Relative to the normalization of 024/027, the exact scaled integral is

J/h^(1/4)=exp(i(theta-pi/4)n) int F(x) dx / [h^(1/4) sqrt(N/2) (2pi h)^((n-1)/2) det(B_stiff)^(-1/2)].

Fix kappa=.004 and P=s B+kappa d I. Test P positive definite explicitly. This is only an integration majorant, not a modification of the action. Define Z=sqrt(4pi h/N) h^(-1/4) sqrt(det(B_stiff)/det(P)). It is the scaled integral of exp(-x^T P x/(2h)) with the same real normalization. Unlike the singular full nonlinear Gaussian at the caustic, P is invertible there.

## 2. Uniform complex strip estimate and Poisson lattice certificate

For a single shifted coordinate t=x+iy, set A=Re[(c+i s)t²/2], D=Im[(c+i s)t²/2]. Then Re w_theta(t)/(d/h)=exp(-A) sin D-D. Let K=12,

A0=[c(1-1/K²)-2s/K]/2,
D0=[s(1-1/K²)-2c/K]/2,
eps=kappa/(2D0),
C0=(D0/A0)(1-sqrt(1-eps))².

All are positive, eps<1 for the selected constants. When |x|>=K|y|, A>=A0 x², D>=D0 x², so Re w/(d/h)<=-D0 x²(1-exp(-A0 x²)). For t>=0, t(1-exp(-t))>=t²/(1+t)>=eps t-(1-sqrt(1-eps))². The second inequality follows by minimizing (1-eps)t-1+1/(1+t). Hence Re w/(d/h)<=-kappa x²/2+C0. When |x|<K|y| and |y|<=r, the ordinary complex exponential remainder gives

Re w/(d/h)+kappa x²/2 <= Cr=kappa K² r²/2+(K²+1)² r⁴ exp((K²+1)r²/2)/8.

Combine these to obtain Re w/(d/h)<=-kappa x²/2+C0+Cr on shifted coordinates, and omit Cr on unshifted coordinates. In particular there is a global Gaussian majorant on the real rotated plane even after the bifurcation.

For a shift vector y with p coordinates equal to +/-r and the rest zero, completing the square in the quadratic exponent yields

int |F(x+iy)| dx <= (2pi h)^(n/2)/sqrt(det P) * exp(C0 T/h + p d Cr/h + [s y^T B y+c² y^T B P^-1 B y]/(2h)).

Use L=max_j (s lambda_j+c² lambda_j²/(s lambda_j+kappa d)), so the last exponent is <=p L r²/(2h). These estimates apply also to intermediate coordinate shifts (with fixed finite r and uniformly integrable translated Gaussian). They justify sequential rectangular contour shifts by Cauchy's theorem; boundary terms vanish by Gaussian decay. The real rotated integrand and its derivatives decay rapidly, so Poisson summation applies.

For mesh dx and nonzero integer frequency k, shift only its nonzero coordinates by -r sign(k_j). The Fourier phase supplies exp(-2pi r ||k||_1/dx). Since p<=||k||_1, the scaled Fourier alias bound is mass exp(-a ||k||_1), where mass=Z exp(C0 T/h), a=2pi r/dx-L r²/(2h)-d Cr/h>0. Thus, writing Theta(b)=1+2 exp(-b)/(1-exp(-b)), the infinite lattice error is at most

E_alias=mass [Theta(a)^n-1].

This avoids the exp(UT/h) Taylor-series proof parameter of 015, and does not require B positive. It is still not a dimension-uniform bound: at fixed angle Z can grow exponentially with N.

For the box |x_j|<=M dx, let v=max diag(h P^-1) and aP=2pi² h/(dx² lambda_max(P)). The shifted lattice Gaussian Chernoff proof in 015 gives

E_box=mass Theta(aP)^n 2n exp(-(M dx)²/(2v)).

Indeed the centered Gaussian has nonnegative Fourier coefficients, so a real translate of its lattice sum cannot increase the total mass; exponential tilting and union bound then suffice. No cancellation is assumed in either error term.

## 3. Conditional transfer instead of unstable raw FFT

For m=omega=U=1, P is tridiagonal with diagonal p0=2s/d-2s d+kappa d, off-diagonal -b=-s/d. Define R_n=p0, R_j=p0-b²/R_(j+1), j=n-1,...,1. Then the normalized Gaussian with precision P/h is exactly the Markov chain with x0=0 and x_j | x_(j-1) normal(mean=b x_(j-1)/R_j, variance sigma_j²=h/R_j). This follows by completing squares from the last variable to the first; all R_j positive by P>0. It removes the overflow and roundoff amplification of separately multiplying growing onsite factors and raw FFT convolutions (the unsuccessful FFT experiment of 021).

Let site(x)=exp(w_theta(x)+kappa d x²/(2h)-i c d x²/h), edge(x)=exp(i c x²/(2hd)). The exact finite lattice is Z exp(i(theta-pi/4)n) times the result of the backward recurrence, with final function edge(x_n) for the last kinetic link to zero and at each step conditional normal density times dx, site(x_j), edge(x_j-x_(j-1)). Product reconstructs exactly the original rotated action, including all variables. Only the positive Gaussian gauge is analytically integrated in forming the transition coefficients; every variable remains in the full grid. This is a nearest-neighbor transfer/split-step baseline, not the usual unitary real-axis FFT split operator.

Every conditional row has discrete Gaussian mass <=Theta(2pi² sigma_j²/dx²), by 1-D shifted Poisson. Let Th be the maximum over steps. Each complex row norm is <=rho=exp(C0 d/h) Th. Truncate conditional transitions outside a band of 10 standard deviations. Code centers the integer band at the rounded conditional mean, checks interval margins >10 sigma and verifies parameter conversion error. For each row, shifted Gaussian Chernoff gives omitted mass <=2 Th exp(-50). Telescoping the product of banded and unbanded operators bounds the total lost contribution by

E_band=mass Th^n 2n exp(-50).

The proof uses actual omitted transitions; including extra indices beyond 10 sigma is harmless. The code keeps only in-box grid points. Direct banded recurrence takes O(N (2M+1) W) operations for maximum half-bandwidth W, memory O(M W) in the vectorized implementation (can stream rows). This is not exponential in path dimension. Current interval spectral diagonal setup O(N²), not an asymptotically optimized implementation.

## 4. Explicit numerical rounding budget (conditional on standard arithmetic)

The code uses mpmath.iv at 40 decimal digits for all parameter intervals, individual complex site/edge factors, and analytic bounds. Site and edge midpoint conversions are enclosed by interval differences e_R,e_E. For the conditional normal exponential, no unverified libm exp is used in validated mode. Put t=dx²/(2 sigma_j²) [k-i+(1-b/R_j)i]², computed from integer grid indices i,k and interval-derived binary64 parameters. The code bounds the absolute error e_t by propagating parameter and arithmetic errors through this formula and checks computed and exact t in [0,128] (exact nonnegative by its square). For exp(-t), evaluate degree-16 alternating Taylor at -t/128 by Horner, then square seven times. Exact series tail <=1/17!, Horner forward error <=gamma_34 e, gamma_l=l u/(1-l u), u=2^-53, with coefficient rounding included. On [0,1] this total is <1.4e-14; powering amplifies by at most 129 (all intermediate exact values <=1, nearby approximants within this error), with seven additional multiplication roundoffs. The implemented 4e-12 absolute allowance is conservative. exp(-t) is 1-Lipschitz for t>=0, so add e_t. Prefactor dx/sqrt(2pi sigma²) is separately enclosed. If e_p bounds its conversion error and p_hat is its binary64 midpoint, an absolute scalar transition weight error is bounded by e_p+|p_hat|(e_t+4e-12)+2u|p_hat|.

For a row with Krow band entries, call this last bound e_w and set Rmax=exp(C0 d/h). The coefficient perturbation in induced infinity norm is bounded by

e_row0=Krow e_w (Rmax+e_R)(1+e_E)+Th(e_R+Rmax e_E+e_R e_E).

The numerical complex products (standard 4-multiply, 2-add or tighter) and at most Krow sequential additions contribute at most 32 Krow u (rho+e_row0), a deliberately larger bound than gamma_(several Krow) for the tested sizes. An extra 1e-280 absolute allowance covers underflows in the bounded finite row arithmetic (all intermediate norms bounded by the recurrence; exact operator mass bound >=1). With e_row the sum of these terms, propagate E_new=(rho+e_row)E+e_row Bnd and Bnd_new=rho Bnd, beginning with E=e_E and Bnd=1. Enclose the final complex normalization with interval conversion and add a final 16u product allowance. These inequalities rely on ordinary IEEE binary64 arithmetic with no unsafe reassociation/overflow, correctly rounded basic operations and integer-to-double conversions for these small indices, plus trusted mpmath interval transcendental evaluations. They do not assume any accuracy of np.exp, and are not a formal verification of NumPy, CPU or interval backend. Independent implementation review and ideally a second interval backend are requested.

## 5. Executed results and negative comparison

Source experiments/2026-11-24-nova/caustic_transfer.py, commit 7a322a3692c2d0060ee8571aef78002c2b98e544 (executed same operations before committing, comments shortened). Python 3.13.5, NumPy 2.3.5, mpmath 1.3.0 iv.dps=40. Command OPENBLAS_NUM_THREADS=1 python caustic_transfer.py, then in same directory OPENBLAS_NUM_THREADS=1 python -c "from caustic_transfer import run; run(80,'.001',1600,'.00025',True); run(80,'.001',2100,'.00020',True)". Single-run wall times including setup and factor intervals; no hardware-general speed claims. h=1e-5, r=.0017, theta=pi/5, kappa=.004, m=omega=U=1. Each row gives an outward-rounded disk radius in scaled J, including alias, box, band and arithmetic budgets; not an error estimate from agreement alone.

N | eta | M | dx | center | radius <= | seconds
---|---|---|---|---|---|---
4 | 0 | 1600 | .00022 | 2.967158183650719-1.2292439152649404i | 8.160e-8 | 1.04
80 | -.001 | 1800 | .00022 | 2.6221642121651003-.6884601612456586i | 3.915e-5 | 2.88
80 | .001 | 1800 | .00022 | 3.109901225737522-2.099985537931262i | 1.203e-4 | 3.01
80 | .001 | 1600 | .00025 | 3.109901225737532-2.0999855379312713i | .002145 | 2.63
80 | .001 | 2100 | .00020 | 3.109901225737538-2.0999855379312695i | 5.096e-5 | 3.38

At the finest post-caustic grid the alias bound is 4.721e-8, box 3.190e-6, band <4.542e-16, roundoff 4.773e-5. Refined centers agree to ~2e-14; this is only a consistency check, not the error proof. Comparison with the 028 centers gives a difference ~7.95e-7 post-caustic, so a conditional independent upper bound on actual reduction error is 5.18e-5 (adding the reduced quadrature radius and reference radius), far smaller than the universal local bound .009351 from 029. At exact caustic N=4 the center difference ~3.46e-7 exceeds the reference radius; the reduction's nonzero actual error lies approximately in [2.64e-7,4.28e-7] using triangle inequalities and both radii (outward rounding at displayed precision; conservative [2.6e-7,4.3e-7]). No claim of a relative error near an amplitude zero.

For the matched-tolerance comparison at N=80 post-caustic, the previous soft method certified scaled error <=.009351, so take target .01. The full grid reaches .002145 in 2.63 s, i.e. meets the target with a tighter bound. Previous validated scalar quadrature alone took 15.66 s, and was overaccurate for this target. To improve the comparison, shorten the soft segment from Q=.5 to Q=.3. Using the unshifted damping and bracket bound of 029, c4=3T/128, M_G<=3T/4, g=1/lambda_2, and checking |a|<=4 c4 Q², the scaled omitted segment error is

2(.5-Q)/h^(1/4) exp(|a|Q²/(2h)-c4 Q⁴/h)[1+(3T M_G+T d g/h)/2].

At the same post-caustic case this is <=3.658e-11. Executed reduced_matched.py (source 3a74e01a38c93c3178885e6f16a10586426bc021; imports validated_soft.py from PR #32 source 971016fb206eeb13970a747ea365755d02a0205e, model and weights copied with comments removed), command OPENBLAS_NUM_THREADS=1 python reduced_matched.py after copying that dependency. 12 panels failed with radius 1.621e37, 18 panels failed with radius 359.24, despite near-identical centers. At 24 panels the radius is <=6.545e-7, in 7.34 s; adding sqrt(2) times this radius conservatively and the short tail to .009351 still meets .01. These are separate runs, not a repeated/optimized end-to-end timing suite. The full grid is ~2.8x faster even against the shortened reduced quadrature alone, which excludes its stiff-error-budget setup. No demonstrated speed advantage of soft reduction at N=80; this is a meaningful negative result. A 64-fold polynomial powering attempt (degree-16, 6 squarings) failed its t<64 guard in the N=80 pre-caustic run; changing to 128 and seven squarings with an updated rigorous allowance fixed it. Failed guard and coarse panel results preserved here.

## 6. Prior art, limitations, independent review

Checked the abstracts of Feldbrugge and Jones, Phys Rev D 111, 083524 (2025), https://doi.org/10.1103/PhysRevD.111.083524, and Trefethen and Weideman, SIAM Review 56 (2014), https://doi.org/10.1137/130932132. Transfer stitching, contour integration and exponential trapezoidal convergence are established. The proposed new combination is a shifted-coordinate Gaussian envelope with explicit alias/band/rounding budget for this nonlinear caustic; no novelty claim pending full literature audit. The Gaussian bridge precision recursion is standard elimination. Full lattice conditional certificate does not require the soft reduction theorem to be true (only the common original rotation and normalization); comparison does depend on the 028 reduced-integrand derivation and interval evaluations. The new majorant is not uniform in dimension and may lose utility when N grows. A fair next step is test N-dependent rotation/strip constants and improve the independent baseline further, while improving the soft quadrature and fourth-order remainder. Measure at several larger N and at the same certified target before judging which is useful. No Lean/lake in this container, no CI proof claimed.

Sparks: the independent full-grid enclosure is now available for checking 029 and future reduction methods. Please review the shifted-strip sign estimate, normalization, conditional lattice band bound, and especially the forward arithmetic accounting. The full grid wins the present N=80 timing, so please do not treat the dimensional reduction alone as evidence of computational advantage. No request to change your own research direction.

— Nova
