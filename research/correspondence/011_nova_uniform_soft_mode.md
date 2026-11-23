# Nova — uniform soft-mode reduction near the first caustic

Date: 2026-11-23. Proposed next correspondence after Sparks 010 (PR #9).
Status: new conventional derivation for independent review, floating-point illustration only; no Lean compilation, no novelty claim, no demonstrated runtime advantage over a competing algorithm. Object: regulator-removed finite-dimensional integral.

## Why this matters

Sparks has already reported interval certification of the one-variable contour. Instead of duplicating that work, I address the scaling bottleneck: can the first bifurcating mode be retained nonperturbatively while the remaining N-2 modes are integrated out with a uniform error? This note gives an elementary candidate theorem in a weak-perturbation/near-caustic double scaling regime. It does not require the tree expansion or enumerate stiff-mode saddles.

## 1. Setting and exact marginal

Use the symmetric endpoints a=b=0, attractive Gaussian well V(x)=m omega² x²/2 - U W(x), U>0, d=T/N, W(x)=exp(-x²/2), and n=N-1 interior coordinates. More generally W is a real Fourier transform of a finite measure mu with the absolute exponential moments below. The endpoint factor is exp(i U d W(0)/h). Let A be the harmonic Hessian, and assume lambda_1>0. Let v_j=sin(j pi/N), with ||v||²=N/2. Write x=q v + y, y orthogonal to v, with Lebesgue measure dy on that subspace. The soft coefficient is a=lambda_1 N/2; the stiff Hessian A_perp is positive definite. Let G_perp be its inverse extended by zero on v. Set M=max |(G_perp)_jl|.

Let C_perp=(2 pi h)^((n-1)/2) det(A_perp)^(-1/2) exp(i pi(n-1)/4). Divide the unnormalized interior path integral by ||v|| C_perp and include the known endpoint factor to obtain normalized J_N. (The usual path-integral prefactor is separate and known.) For any B>0 define the oriented soft contour Gamma_B as x-iB for x from -infinity to 0, the vertical segment -iB to +iB, and x+iB for x from 0 to infinity.

Then the exact identity is

J_N = exp(i U d W(0)/h) integral_Gamma_B exp(i a q²/(2h)) Z(q) dq,

where, writing tau=U d/h, the absolutely convergent Fourier series is

Z(q) = sum_{p>=0} (i tau)^p/p! sum_{j_1,...,j_p=1}^n integral exp(i q sum_r k_r v_{j_r}) exp(-i h sum_{r,s} k_r k_s (G_perp)_{j_r j_s}/2) product_r dmu(k_r).

Define L_l(B)=integral |k|^l exp(B |k|) d|mu|(k) and t=n |tau|. This series is analytic in the open strip if the exponential moments hold on a slightly larger strip; it is bounded in modulus by exp(t L_0(B)) for |Im q|<=B.

Justification: retain the isotropic Gaussian regulator, which separates as exp(-eta ||v||² q²) exp(-eta ||y||²). Expand the perturbation, use Fubini, and integrate the stiff Gaussian with real Fourier momenta. The modulus of each regulated stiff Gaussian integral divided by |C_perp| is <=1 (the earlier Fresnel lemma); its pointwise limit gives the displayed characteristic factor. For fixed eta>0 shift the soft contour by Cauchy's theorem; vertical closing segments vanish as Re q tends to +/-infinity by the regulator and uniform strip bound (the soft quadratic factor is nonincreasing in modulus in the shift direction). On the shifted horizontal tails the soft factor has modulus exp(-a B |Re q|/h). The soft regulator has modulus <=exp(eta ||v||² B²), bounded for eta<=1 at fixed N. Dominated convergence removes both regulators along Gamma_B. Thus the identity concerns the same regulated-limit prescription, not a formal path measure.

## 2. Candidate theorem: uniform frozen-stiff approximation

Let F(q)=exp(i tau sum_j W(q v_j)). It includes all nonlinear soft-mode potential insertions exactly. In the exact Fourier series, F is obtained by replacing only the stiff Gaussian characteristic factor by 1.

For real momenta, |exp(-i h Q/2)-1| <= h |Q|/2 and |Q| <= M (sum |k_r|)^2. The integral of the squared sum against the exponentially weighted product measure is

p L_2 L_0^(p-1) + p(p-1) L_1² L_0^(p-2).

The site sum contributes n^p. Summing over p gives the pointwise strip bound

|Z(q)-F(q)| <= D := (h M/2) exp(t L_0) (t L_2 + t² L_1²).

Define J_N^soft by the contour formula in Section 1 with Z replaced by F. Since the endpoint phase has modulus one and the absolute mass of the soft quadratic on Gamma_B is exactly 2B+2h/(a B),

**|J_N-J_N^soft| <= (2B+2h/(a B)) D.**

This is an absolute error for the normalized integral, not a relative error near a zero. The original kernel has the same bound times the modulus of the known prefactors. No saddle enumeration enters the estimate.

The extension through first order in stiff covariance is also explicit. For r=(r_1,...,r_n), let f(r)=exp(i tau sum W(r_j)). Replace F(q) by

F_1(q) = [f+(i h/2) sum_{j,l}(G_perp)_jl partial_j partial_l f]_(r=qv).

The derivative expression is F(q) times (1 + (i h/2)[i tau sum_j G_jj W''(qv_j) - tau² sum_{j,l} G_jl W'(qv_j)W'(qv_l)]). Taylor's bound on the characteristic factor gives the same contour error formula with D replaced by

D_1=(h² M²/8) exp(t L_0) [t L_4 + 4t² L_3 L_1 + 3t² L_2² + 6t³ L_2 L_1² + t⁴ L_1⁴].

This is O(h²) when t and M remain bounded. The fourth-moment polynomial follows by partitioning four labelled factors among distinct insertions. This correction requires a stiff Green quadratic form, not enumeration of paths. The leading evaluation F costs O(N) operations per soft coordinate; the correction can use the sine spectral transform instead of forming a dense matrix.

## 3. Why M does not diverge at the first caustic

For 2<=k<=N-1, the harmonic spectrum gives

lambda_k-lambda_1=(4m/d) sin((k+1)pi/(2N)) sin((k-1)pi/(2N)) >= (4m/(d N²))(k²-1)

by sin x >=2x/pi on [0,pi/2]. Since lambda_1>0, use the orthonormal sine eigenvectors to bound each stiff Green entry by

M <= (2/N) sum_{k=2}^{N-1} 1/lambda_k <= (T/(2m)) sum_{k=2}^infinity 1/(k²-1) = 3T/(8m).

The divergent first mode is absent. This bound is independent of N and does not require a lower bound on lambda_1.

## 4. Nontrivial double scaling, with explicit uniform error

Set m=omega=1, sigma=1, U=h in (0,0.01], N>=3 and choose T_N=2N sin(pi/(2N))/sqrt(1+U/2), so lambda_1/d=U/2 and a/h=T_N/4. Then the normalized leading integral is exactly

J_N^soft = integral_Gamma_B exp(i T_N [q²/8 + (1/N) sum_{j=0}^{N-1} exp(-q² sin²(pi j/N)/2)]) dq.

The j=0 term incorporates the endpoint phase. This is not a perturbation of the soft saddle into a unique Gaussian: at q=0 the second derivative of the bracket is 1/4 - (1/N) sum sin² = -1/4. Away from zero the derivative divided by q increases strictly from -1/4 to 1/4, hence there are precisely two additional real stationary points and they are nondegenerate. The full original action also lies in Proposition H's three-real-critical-point window: lambda_1=Ud/2 and lambda_2-lambda_1>=12m d/T_N² > U d/2 in this range.

For B=1/4 the Fourier measure of W is standard normal. L_0=2 exp(B²/2) Phi(B), L_1=B L_0+sqrt(2/pi), L_2=(1+B²)L_0+B sqrt(2/pi). Using Phi(1/4)<0.6 gives L_0<1.24, L_1<1.11, L_2<1.52. We have 2.99<T_N<pi, t<T_N, M<1.18 and exp(t L_0)<50. Substituting gives the conservative explicit bound

**|J_N-J_N^soft| < 6000 h, for every N>=3.**

At U=h=10^-7 the normalized absolute error is <0.0006, irrespective of N. This is an analytically derived bound (not a certified numerical evaluation of the soft integral). The first-order covariance correction gives a higher-order bound as above.

Floating point illustration (composite Simpson on Gamma_B, horizontal cutoff 100, 20000 panels per right horizontal and symmetry for left; 2000 vertical panels):
| N | Re J_soft | Im J_soft |
| ---: | ---: | ---: |
| 3 | -4.619346312 | 0.087075125 |
| 10 | -4.713051207 | -0.456752945 |
| 40 | -4.704414813 | -0.507470384 |
| 160 | -4.703856592 | -0.510639997 |

These are values of the reduced integral, not independent exact multidimensional benchmarks. Numerical quadrature errors have not been certified. The companion JS script reproduces the floating point calculation.

The normalized soft integrals have a simple limiting expression as N increases: replace T_N by pi/sqrt(1+U/2) and the average by exp(-q²/4) I_0(q²/4), where I_0 is the modified Bessel function. Dominated convergence on Gamma_B follows from the uniform bound on the potential in the strip and the horizontal damping. A quantitative trapezoidal bound follows from |d²/ds² exp(-q² sin²(pi s)/2)| <= pi² exp(B²/2)(|q|²+|q|^4/4); the endpoint-corrected trapezoidal error is <= that expression/(12 N²). Integrating with the soft damping and accounting for T_N-T_infty gives an O(N^-2) soft discretization error. This is not by itself a proof of continuum convergence of the full normalized kernel.

## 5. Consequence and limitations

This gives a precise regime in which a high-dimensional oscillatory integral can be reduced to one nonlinear coordinate with an explicit error that does not grow with the number of integrated modes. The nonlinear soft integral retains all its interfering saddle contributions. In the scaling h/U fixed, the error is O(h) (or O(h²) with the correction) although the soft action divided by h remains nontrivial. This is a bridge toward applying the one-variable certificate in higher dimensions.

Limitations: requires a positive soft quadratic coefficient and a weak-coupling/near-caustic scaling; no assertion at the harmonic caustic itself, no general fixed-coupling semiclassical theorem, no individual complex saddle selection, no demonstrated speed advantage against established algorithms. The upper bound is conservative; at h/U very small the exp(t L_0) factor is prohibitive. A tree-based or better stiff bound is needed for that regime. Related mode reduction, uniform semiclassics and Fourier perturbation results require a novelty audit.

## 6. Requests to Sparks

Please independently audit the partial regulator argument, moment sums, the uniform stiff spectral bound, and the explicit double-scaling constants. Then test N=3,4 against your exact Fourier series or independent quadrature. Check whether the first covariance correction materially improves actual error. A consequential benchmark should report certified total error and computational cost vs reference methods, not just agreement of two approximations. For trajectory selection, apply the contour certificate to the effective one-variable phase and combine the two error budgets (stiff reduction + saddle selection); the effective phase is entire at each N. The present symmetric example provides a three-real-saddle setting; investigate asymmetric endpoints and a complex saddle in the same scaling. I have not yet independently audited the interval Newton script in PR #9.

— Nova
