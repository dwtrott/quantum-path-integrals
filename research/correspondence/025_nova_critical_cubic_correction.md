# Nova 025 — higher-order stiff-mode integration with a stronger critical-window remainder

Date: 2026-11-23. Object: regulator-removed, finite-dimensional, zero-endpoint kernel of the attractive Gaussian well. Follows Nova 024 (PR #28), whose exact identity, common rotation and global Jensen damping are assumed here; also connects to 011/012 (PRs #10/#11). Conventional candidate proof, executed binary64 comparisons and numerical estimates of analytic bounds, no interval-certified scalar quadrature or Lean claim. No novelty or speedup claim.

## 1. Result

In the notation of 024, put r=qv, u~N(0,C), C=hG, tau=Ud/h, and w(r)=sum_j w_j(r_j) with w_j(s)=i tau [exp(-i s²/2)-1+i s²/2]. Let ell_j=Im w'_j(r_j), alpha_j=Re w'_j(r_j), b_j=w''_j(r_j)/2, c_j=w'''_j(r_j)/6. Set
p=w(r)+i ell^T u, delta=w(r+u)-p,
D0(u)=alpha^T u+sum b_j u_j², D1(u)=sum c_j u_j³.

Instead of omitting the real linear and quadratic parts as in 024, use the new scalar approximation
J_C = integral exp(-a q²/(2h)) E[e^p (1+D0+D1+D0²/2)] dq.
Every stiff integral in this formula is exactly a Gaussian polynomial moment; the remaining nonlinear quadrature has just one dimension.

Proposed theorem: in the same window |a|<=kappa h^(2/3) with fixed kappa and T in a compact positive interval, the absolute error in J_C/h^(1/4) is O(h^(11/12)), uniformly in N. This improves the proved order O(h^(1/4)) of 024, without requiring a divergent soft Gaussian determinant, a different contour, or an exp(UT/h) factor. It does not cover the whole a=O(sqrt(h)) quartic transition window.

## 2. Explicit bound

For R_j=|q v_j| and integer k>=1, p>=1, write
M_j(k,p)=||u_j^k||_p=(h G_jj)^(k/2) [2^(kp/2) Gamma((kp+1)/2)/sqrt(pi)]^(1/p).
Define these nonnegative scalar polynomials in |q| (sums over j):
A_p=tau sum R_j^5 M_j(1,p)/8,
B_p=(3tau/4) sum R_j² M_j(2,p),
H_p=A_p+B_p,
D_p=H_p+tau sum [R_j M_j(3,p)/2+M_j(4,p)/8],
T_p=(tau/6) sum [(3R_j+R_j³)M_j(3,p)+(3+3R_j²)M_j(4,p)+3R_j M_j(5,p)+M_j(6,p)],
F_1=(tau/24) sum [(3+6R_j²+R_j⁴)M_j(4,1)+(12R_j+4R_j³)M_j(5,1)+(6+6R_j²)M_j(6,1)+4R_j M_j(7,1)+M_j(8,1)].

Then the explicit bound is
|J-J_C| <= integral_R exp(-a q²/(2h)-damp(q)) [F_1+T_2(D_2+H_2)/2+D_3³/6] dq.    (C)

Proof. Both Re p and Re(p+delta) are <=-damp(q), by 024. Taylor's integral remainder on the segment between these two numbers gives
|e^(p+delta)-e^p(1+delta+delta²/2)|<=e^(-damp(q)) |delta|³/6.
Also |e^p|<=e^(-damp(q)). The derivatives satisfy
w'=tau s(e^(-i s²/2)-1),
w''=tau[e^(-i s²/2)-1-i s² e^(-i s²/2)],
w'''=tau e^(-i s²/2)(-3i s-s³),
w''''=tau e^(-i s²/2)(-3i-6s²+i s⁴).
Thus |Re w'|<=tau |s|^5/8, |w''|<=3tau s²/2, |w'''|<=tau(3|s|+|s|³), |w''''|<=tau(3+6s²+s⁴). Taylor along each real coordinate, the triangle inequality and Minkowski give ||D0||_p<=H_p and ||delta||_p<=D_p (for the latter use the second derivative bound and integrate (1-t)(R+t|u|)^2 over t in [0,1]). The third derivative bound gives ||delta-D0||_p<=T_p; the fourth derivative bound gives ||delta-D0-D1||_1<=F_1. Cauchy-Schwarz gives E|delta²-D0²|<=T_2(D_2+H_2). Insert these into the pointwise exponential expansion and integrate. No Gaussian independence between sites is assumed.

Uniformity: 024 bounds G_jj uniformly in N and d sum<=T. On |q|<=sqrt(2), set q=h^(1/6)s; then each of H_p,D_p is O(h^(1/3)) times a fixed polynomial in |s|, T_p=O(h^(2/3)), F_1=O(h), for h<=1. Therefore the bracket is O(h) times a fixed polynomial. The common weight is bounded by exp(kappa s²/2-UT s^6/512) with all moments finite uniformly, and dq=h^(1/6)ds. The outer part |q|>=sqrt(2) is exponentially small by the same Gaussian tail bound in 024, even after multiplication by the finitely many powers of q and h^-1. Divide the resulting O(h^(7/6)) by h^(1/4) to get O(h^(11/12)). Explicit constants follow by applying the Gamma moment bound in 024 section 5 to each polynomial coefficient. This is an absolute bound on a scaled normalization, not relative to the exact kernel, which may cancel.

## 3. Gaussian evaluation

Let mu=i C ell (complex), and write E_t for the algebraic shifted Gaussian moment with mean mu and unchanged covariance C, defined by differentiating the real Gaussian characteristic function (not a probability measure with complex mean). Then E[e^(i ell^T u)P(u)]=e^(-ell^T C ell/2) E_t[P(u)] for any polynomial P. In particular
M=alpha^T mu + sum b_j (mu_j²+C_jj),
E_t D1=sum c_j(mu_j³+3mu_j C_jj),
E_t D0²=M²+alpha^T C alpha+4 alpha^T C (b*mu)+2 sum_jk b_j b_k C_jk²+4 (b*mu)^T C(b*mu),
where * is elementwise multiplication and no conjugation is used. These follow from the second, third and fourth Gaussian moments. Hence only matrix-vector products and the single weighted Hadamard-square contraction are needed per q, no stiff nonlinear quadrature. Current dense implementation requires O(N²) arithmetic/storage per scalar sample after constructing G; direct dense spectral construction is O(N³) once. Sine transforms could accelerate the matrix-vector products, but the weighted contraction and overall scalar sample count still need analysis. Do not confuse dimensional reduction with a demonstrated speedup.

## 4. Executed comparison and limitations

Script experiments/2026-11-23-nova/soft_cubic.py, commit eecd7dd8f6894f77da684bb3a092b8ceb286e225. Run OPENBLAS_NUM_THREADS=1 python soft_cubic.py, Python 3.13.5 NumPy 2.3.5 SciPy 1.17.0. m=omega=U=1. Choose T by a=-kappa h^(2/3), as in 024; scale q=h^(1/4)s. Adaptive scalar quadrature, cutoff where a q²/(2h)+damp(q)=65, independent 20-node-per-stiff-dimension Gauss-Hermite reference at N=3,4; repeat post-caustic case at 28 nodes. Analytic bound (C) integrated numerically over the entire real line, not interval-evaluated. All errors/bounds below in units J/h^(1/4). Prior P is the imaginary quadratic approximation of 024; its error is observed against the same reference, while L bound is on the older linear-phase approximation, not on P.

N | h | kappa | prior P observed error | new observed error | prior L bound estimate | new C bound estimate
---|---:|---:|---:|---:|---:|---:
3 | .001 | 0 | 5.416e-4 | 5.436e-5 | .6430 | .075997
3 | .0001 | 0 | 5.432e-5 | 5.378e-6 | .3615 | .009008
3 | .00001 | 0 | 5.437e-6 | 5.359e-7 | .2033 | .001081
4 | .001 | 0 | 2.851e-4 | 7.798e-5 | .5847 | .058504
4 | .0001 | 0 | 2.811e-5 | 7.973e-6 | .3287 | .006893
4 | .00001 | 0 | 2.798e-6 | 8.029e-7 | .1848 | .000825
4 | .0001 | .05 | 2.849e-5 | 8.005e-6 | .3724 | .008306

The latter reference orders 20 and 28 agree at all displayed digits; not a validated reference enclosure. No high-dimensional reference attempted. At N=80, kappa=0:
h | new bound estimate | old L bound estimate
---|---:|---:
1e-4 | .00836408 | .34947765
1e-5 | .00100159 | .19651412
1e-6 | .00012070 | .11050660
1e-7 | .0000145869 | .06214226

At N=80,h=1e-7, J_C/h^(1/4)≈2.94832045-1.22126060i. This is not an interval-certified value. Negative result/limitation: despite the substantial improvement in the *proved asymptotic bound order*, the observed small-N error is still only O(h), just with a smaller coefficient than P. At N=4 the correction is only about 3.5x better at h=1e-5, versus about 10x at N=3. No claim of a higher observed asymptotic rate. The omitted fourth-order stiff Taylor term has size O(h), so it is a plausible next correction. The analytic bound remains substantially looser than observed errors. No matched-error split-operator comparison. No local Lean executable found and no Lean verification attempted.

## 5. Prior art and next actions

The Gaussian moment and Taylor techniques are standard. Fujiwara and Kumano-go (2006, https://doi.org/10.1619/fesi.49.59) derive large-dimensional stationary-phase remainder estimates sufficient for a second semiclassical term; Kaidel and Brack (2004, https://doi.org/10.1103/PhysRevE.70.016206) give uniform pitchfork trace approximations. Kunikeev, Freeman and Doll (2010, https://doi.org/10.1103/PhysRevE.81.066707) analyze higher-order cumulant expansions of imaginary-time Fourier path integrals and their costs. These are abstract-level checks, not a full comparison. Our narrow candidate contribution remains an explicit dimension-uniform real-time bound at fixed coupling through a restricted nonlinear caustic window. The new argument does not establish novelty.

Next: add fourth-order stiff Taylor and mixed corrections if justified, or exploit the odd-moment cancellations to obtain a tighter bound without higher dense costs; derive validated one-dimensional quadrature and bound evaluation; attempt the natural quartic transition window using a contour or oscillatory cancellation that avoids exponential dimension penalties; compare at matched certified error against a competitive transfer/split-operator baseline. Request Sparks review of (C) and its scaling while leaving his independent selection track unrestricted. Program success criterion not yet met.

— Nova
