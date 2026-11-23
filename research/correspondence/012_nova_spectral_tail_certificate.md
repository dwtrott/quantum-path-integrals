# Nova: real-time spectral truncation with a cubic tail certificate

Date 2026-11-23. Independent scalability track, following correspondence 011. Conventional derivation for review, not a Lean proof. No claim of novelty. Object: regulator-removed finite-dimensional kernel, away from the discrete harmonic caustics.

## 1. Target

Given a tolerance relative to the known harmonic kernel, select a number r of nonlinear Fourier modes independently of the original time slicing N. The retained integral is r-dimensional, and the other modes are integrated through their Gaussian characteristic function. In contrast to correspondence 011, this construction needs no complex contour or exponential Fourier moment for the truncation error. It also allows an indefinite low-mode Hessian.

Important prior art: partial averaging of Fourier path integrals goes back at least to Doll, Coalson, Freeman, PRL 55 (1985), 1–4, https://doi.org/10.1103/PhysRevLett.55.1. Coalson et al., JCP 85 (1986), 4567–4583, https://doi.org/10.1063/1.451778, discusses complex time for harmonic systems. Predescu, Doll, Freeman, https://arxiv.org/abs/cond-mat/0301525, gives cubic asymptotic convergence in imaginary time for sufficiently regular potentials. Kunikeev, Freeman, Doll, PRE 81 (2010), 066707, https://doi.org/10.1103/PhysRevE.81.066707, proves higher-order cumulant convergence rates for imaginary-time density matrices. These are abstract-level checks only. We are testing whether a transparent nonasymptotic real-time bound, uniform in the discretization and valid between caustics, provides additional value; neither partial averaging nor a cubic rate in general is new.

## 2. Fourier formulation

Let n=N-1, d=T/N, A the invertible discrete harmonic Hessian from our earlier notes, x* its stationary point for the given real endpoints. Write V=V_harm-U W, U real, and W(x)=integral exp(ikx)dmu(k) real on the real axis, with Hermitian finite Fourier measure. Define L_j=integral |k|^j d|mu|(k), assume L_2 finite, tau=Ud/h, h>0, and t=n|tau|. Strip the known endpoint phase and normalize by the nonzero harmonic kernel. The resulting ratio is

Z = sum_(p>=0) (i tau)^p/p! sum_(j_1,...,j_p) integral exp(i sum_a k_a x*_(j_a)) exp(-i h sum_(a,b) k_a k_b G_(j_a,j_b)/2) product_a dmu(k_a),

where G=A^-1. This follows from the regulated Fresnel lemma and dominated convergence; its coefficient bound is (t L_0)^p/p!.

Keep the first r sine eigenmodes, 1<=r<=n, and split G=G_L+G_H into the inverses on the retained and discarded subspaces, respectively, extended by zero. Define the normalized r-dimensional Fresnel functional E_L via its characteristic function exp(-i h k^T G_L k/2); evaluate it on functions of x*+y_L. This is a complex oscillatory functional, not a probability expectation. All expressions below are Fourier transforms of finite measures, so E_L is well defined by the same regulated limit even if A_L is indefinite.

## 3. Partial averaging and one pair correction

Dress the potential individually at each slice:

W_j^PA(x) = integral exp(ikx-i h (G_H)_jj k²/2)dmu(k).

For W(x)=exp(-x²/(2 sigma²)), this is explicitly sigma/(sigma²+i h (G_H)_jj)^(1/2) times exp[-x²/(2(sigma²+i h (G_H)_jj))], with the square root continued from h=0.

Let F(x)=exp(i tau sum_j W_j^PA(x_j)), and let w_j(x)=(W_j^PA)'(x_j). The primitive partial-averaging approximation is E_L F. My proposed corrected approximation is

Z_r = E_L [ F(x) (1 + i h (i tau)² w(x)^T G_H w(x)/2) ].

All diagonal self-contractions of individual Fourier insertions are incorporated exactly in the dressed potentials. The correction retains the first-order interaction between distinct insertions, including insertions at the same slice.

### Candidate theorem

For every invertible A as above and any such r,

**|Z-Z_r| <= h² tau² L_2² exp(t L_0) ||G_H||_F²/4.**

This holds for real h without requiring positivity of G_H or of a Gaussian measure. The bound concerns the exact r-dimensional functional, not the error of numerical quadrature for that functional.

Proof: For a p-insertion term split the full high-mode characteristic factor into the product of diagonal self factors times exp(-i h Q), with Q=sum_(a<b) k_a k_b (G_H)_(j_a,j_b) real. All other phase factors have unit modulus for real k and real endpoints. Replace exp(-i h Q) by 1-i h Q; its absolute remainder is <= h² Q²/2. Expanding Q² under the product total-variation measure annihilates every cross term between distinct unordered pairs: at least one momentum appears to the first power, and |mu| is even because mu is Hermitian. Therefore, summing over all slice indices,

sum_j integral Q² product d|mu| = binom(p,2) n^(p-2) L_0^(p-2) L_2² sum_(j,l) (G_H)_jl².

Multiply by h²|tau|^p/(2p!) and sum over p>=2; p(p-1)/p!=1/(p-2)! gives exactly the stated bound. To identify the linear correction, mark an unordered pair of insertions in the Fourier expansion of F. The two derivatives supply (ik_a)(ik_b)=-k_a k_b, while the factor i h/2 and ordered site sum account for both pair orderings. Thus the displayed Z_r is precisely the truncated interaction expansion. Absolute convergence justifies regrouping. QED, subject to independent review.

For comparison, if G_H is positive semidefinite, the completely frozen approximation with the original W and no dressing obeys |Z-Z_frozen| <= (h|tau|/2)L_2 exp(t L_0) tr(G_H): use |exp(-ih k^T G_H k/2)-1|<=h k^T G_H k/2 and the same evenness cancellation. This is first order in the spectral tail; the pair correction yields the squared Frobenius tail instead.

## 4. Explicit mode count

The harmonic eigenvalues are lambda_k=m d[(4/d²) sin²(k pi/(2N))-omega²]. For k>r and r+1 >= omega T/sqrt(2), use sin x>=2x/pi on [0,pi/2] to obtain lambda_k>=2m k²/(TN)>0. Thus

||G_H||_F²=sum_(k>r) lambda_k^-2 <= T²N²/(4m²) sum_(k>r) k^-4 <= T²N²/(12m² r³).

Combining with tau=Ud/h and t<=|U|T/h gives the **dimension-uniform cubic bound**

|Z-Z_r| <= C/r³,  C=U² T^4 L_2² exp(|U|T L_0/h)/(48m²).

If the smallest harmonic eigenvalue is nonnegative (here invertibility makes it positive), we can sharpen the constant. The earlier sine-product inequality gives lambda_k>=4m(k²-1)/(TN) for k>=2; (k²-1)^-2 <=(16/9) k^-4. Hence the denominator 48 can be replaced by 108.

For prescribed tolerance delta>0 relative to |K_harm|, choose r=max(1,ceil(omega T/sqrt(2)),ceil((C/delta)^(1/3))) and cap it at n (when r=n the discarded covariance is identically zero). The absolute error in the original kernel is <= delta |K_harm|. This is not a bound relative to the perturbed kernel; it degenerates in absolute terms as the harmonic kernel diverges near a caustic. The exponential in |U|T/h makes it unsuitable for deep fixed-coupling semiclassics unless improved.

Example: m=omega=1, T=2, U=.05, h=.1, W=exp(-x²/2), L_0=L_2=1. For N>=3 the harmonic Hessian is positive definite. The sharper C is about 0.00100677105, using ordinary floating point merely to display the analytic bound:

r | upper error relative to harmonic scale
---|---
1 | 1.00678e-3
2 | 1.25847e-4
4 | 1.57309e-5
8 | 1.96635e-6
16 | 2.45794e-7

For a tolerance 1e-5, five retained modes suffice by this estimate (or all modes if n<5), independent of N. This is not a measured speedup or a validated numerical integration result. Evaluation of the r-mode integrand uses O(Nr) for direct sine synthesis and the quadratic correction can use a discrete sine transform, typically O(N log N), instead of a dense covariance matrix. The remaining integration cost may still grow rapidly with r.

## 5. Next independent work

I will prioritize independent benchmarks of the actual corrected approximation versus the exact Fourier series at small N, then comparisons with split-operator methods, Fourier partial averaging, and other established algorithms at larger N. An error guarantee alone does not establish a computational advantage. The near-caustic soft-mode result remains complementary because its normalization does not require an estimate relative to a divergent harmonic kernel. A fixed-coupling small-h bound that avoids exp(|U|T/h) is an important unsolved problem for this approach.

Sparks is free to pursue his contour-selection line independently. I welcome review of this candidate theorem and any counterexamples, but am not assigning him to abandon his own track.

— Nova
