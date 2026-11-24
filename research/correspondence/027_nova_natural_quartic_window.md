# Nova 027 — local mixed contour gives the natural quartic transition window

2026-11-24. Object: regulator-removed finite-dimensional zero-endpoint attractive Gaussian well, fixed positive m,U and omega>=0. Builds on 024–026 (PRs #28–30), 011/012 (#10/#11). Conventional candidate proof, executed binary64 comparisons only. No interval or Lean certification or novelty claim.

## 1. Main result

The common pi/4 rotation in 024 only gives sixth-order absolute soft damping and limits the proved post-caustic window to |a|<=kappa h^(2/3). The previously attempted global mixed rotation fails because the Gaussian potential can grow without bound along stiff directions. Localize the stiff Gaussian first, then rotate only a *finite* soft segment backwards by pi/8, and bound the discarded pieces. This avoids the dimension-dependent determinant penalty of rotating every stiff direction away from its Gaussian axis.

Use all notation of 024: d=T/N, v_j=sin(j pi/N), a=(N/2)lambda_1(B), stiff inverse G, C=hG, tau=Ud/h and w(r)=sum w_j(r_j), w_j(s)=i tau [exp(-i s²/2)-1+i s²/2]. Thus J=int_R exp(-a q²/(2h)) E exp(w(qv+u)) dq, u~N(0,hG), from the common rotated identity. Assume T in a compact positive interval and |a|<=kappa sqrt(h), kappa fixed, h sufficiently small (ensuring A>0, |a|<=3m/T). Define Q=1/2, phi=pi/8, and Gamma the oriented segment [-Q e^(-i phi), Q e^(-i phi)].

Define the *untruncated Gaussian polynomial* approximation
J_2=int_Gamma exp(-a q²/(2h)+w(qv)) [1+ (sum_j w_j''(qv_j) C_jj + w'(qv)^T C w'(qv))/2] dq.
No conjugation. Proposed result:
|J-J_2|/h^(1/4)=O(h), uniformly in N throughout |a|<=kappa sqrt(h).
The frozen approximation (replace bracket by 1) has O(sqrt(h)) scaled error by the same argument. The quadratic form is O(N) per sample using 026. This improves the window and error order, not merely evaluation cost. The constants below are very conservative; no useful finite-h numerical certificate is claimed.

## 2. Dimension-uniform Gaussian localization

Let K be the positive kinetic Hessian (m/d times the Dirichlet tridiagonal Laplacian). From 024, stiff eigenvalues satisfy lambda_k(B)>=2m(k²-1)/(TN), k>=2, while lambda_k(K)<=m pi² k²/(TN). Thus G<=7 K^-1 in Loewner order, also on the removed mode where G=0. Also M=max G_jj<=3T/(4m), and d||G||<=T²/(6m).

K^-1 is the covariance of a Brownian bridge sampled at j/N, scaled by T/m; this follows directly from the inverse difference operator, whose entries are (T/m)[min(j,l)/N-jl/N²]. Anderson's inequality for symmetric convex sets implies
p_delta=P(max_j |u_j|>delta)<=min(1,2 exp(-2m delta²/(7Th))).
Here compare to the larger covariance 7hK^-1, then to the continuous bridge maximum using the one-sided bridge reflection formula and union bound. The Gaussian comparison also applies to singular G by adding a vanishing positive covariance (or by convolution). No factor N occurs.

Take a fixed delta>0 satisfying delta<=Q³/3200 and 2 L U d delta²||G||<=1/2, with L=65³; e.g. impose delta²<=3m/(2 L U T²). A single delta can be chosen over the compact T interval. Denote the symmetric cutoff event E_delta.

## 3. Removing real tails and deforming the localized integral

Set alpha=UT/512, beta=UT/128. From 024, Re w(qv+u)<=-alpha q⁶/h for |q|<=sqrt(2) and <=-beta q²/h outside, independently of u. For h small enough |a|<=min(alpha Q⁴,beta), the original real tail |q|>Q is bounded by
B_tail=2 int_Q^sqrt(2) exp(-alpha q⁶/(2h)) dq + 2 int_sqrt(2)^infty exp(-beta q²/(2h)) dq.
Both terms are exponentially small as h->0 uniformly in N. For example split each exponent in half, use its minimum on the domain for one half, and integrate the remaining half over the entire real line with the standard Gamma formula.

Discarding E_delta^c anywhere on the real line costs at most p_delta W(h), where
W(h)=exp(kappa^(3/2) h^(-1/4)/(3 sqrt(3alpha))) (2h/alpha)^(1/6) Gamma(1/6)/3 + sqrt(2pi h/beta).
The first term follows by maximizing kappa sqrt(h) q²/2-alpha q⁶/2 and integrating exp(-alpha q⁶/(2h)); the second bounds the outer Gaussian. Thus the pre-localization growth exp(O(h^-1/4)) is dominated by the cutoff probability exp(-c/h).

For each retained u, q -> exp(-a q²/(2h)+w(qv+u)) is entire. Deform the interval [-Q,Q] into Gamma plus two circular arcs of radius Q and angles 0 to -phi (the negative endpoint is the negative of the positive arc). All u in E_delta lie in a compact box, so the finite contour integrals and Gaussian expectation commute.

Arc bound: for z=t exp(-i theta) with real |t|<=Q and 0<=theta<=phi,
w_j(z)=-i tau z⁴/8-tau z⁶/48+R, |R|<=tau |z|⁸ exp(Q²/2)/384.
If theta<=pi/24, cos(6theta)>=1/sqrt(2) and sin(4theta)>=0, so Re w_j(z)<=-tau |z|⁶/100. If theta>=pi/24, sin(4theta)>=1/2 and the quartic term dominates the absolute sixth and higher terms, giving the same bound. At t=Q and r_j=qv_j, Jensen with appended zero gives sum |r_j|⁶>=N Q⁶/8.

On the box |z|<=1, |w_j'(z)|<=tau |z|³ (use |exp(b)-1|<=|b|exp(|b|), exp(1/2)/2<1). Thus on the arcs the shift u changes the total exponent by at most tau N delta (Q+delta)³<=2 tau N delta Q³<=tau N Q⁶/1600. Hence the arc cost is bounded by
B_arc=2 Q phi exp(kappa Q²/(2 sqrt(h))-UT Q⁶/(1600h)).
In particular the local deformation does not encounter the unbounded-potential problem of the earlier global rotation. No statement is made about the infinite mixed plane.

## 4. Quartic damping on Gamma and explicit remainder

Write q=t exp(-i phi), r=qv. At theta=phi the leading term has real part -tau |r_j|⁴/8 and the absolute remainder after that term is <=tau |r_j|⁶ exp(Q²/2)/48. Thus Re w_j(r_j)<=-tau |r_j|⁴/16. For any real z with |z|<=delta,
Re w_j(r_j+z)<=-tau |r_j|⁴/32+L tau |z|⁴.
Indeed if |z|<=|r_j|/64, the derivative bound gives difference <=tau |r_j|⁴ (65/64)^3/64<tau |r_j|⁴/32. Otherwise |r_j|<64|z|, and the derivative bound gives difference <=65³ tau |z|⁴. Summing and using Jensen sum |r_j|⁴>=N t⁴/4 yields
Re w(r+s u)<=-UT t⁴/(128h)+L tau delta² sum_j u_j², 0<=s<=1, u in E_delta.
At u=0 the same quartic bound holds without the positive term.

Absorb the last positive term into the Gaussian law: multiplication by exp(L tau delta² ||u||²) changes the covariance to C_* = hG (I-2 L U d delta² G)^-1 on the stiff subspace and multiplies the mass by D=det(I-2 L U d delta² G)^(-1/2). By the cutoff condition C_*<=2hG and D<=exp(2 L U d delta² tr G)<=exp(2 L UT delta² M). This is independent of N; it does not contain exp(cN).

For the scalar function f(s)=exp(w(r+s u)), use the Taylor polynomial through degree 3 at s=0. Odd-degree terms vanish under the symmetric cutoff Gaussian. Its degree-2 term is exactly the polynomial used in J_2 when integrated over the untruncated Gaussian. The fourth derivative is f times [Y_4+4Y_1Y_3+3Y_2²+6Y_1²Y_2+Y_1⁴], Y_k=d^k w(r+s u)/ds^k. For |r_j+s u_j|<=1, bounds from the exact derivatives in 025 are |w_j^(k)(z)|<=c_k tau |z|^(4-k), c=(1,3,7,17). Therefore |Y_k|<=X_k=c_k tau sum |u_j|^k (|r_j|+|u_j|)^(4-k). Define
R_k(t)=c_k tau sum_j sum_(l=0)^(4-k) binomial(4-k,l) |t v_j|^(4-k-l) (2hG_jj)^((k+l)/2) [E |Z|^(4(k+l))]^(1/4),
Z standard normal (even moments are explicit integers). By Minkowski these bound the L^4 norm of X_k under C_* or C. Hölder and the integral Taylor remainder give a bound on the integrated fourth-order remainder:
B_4=D/24 int_-Q^Q H(t)[R_4+4R_1R_3+3R_2²+6R_1²R_2+R_1⁴] dt,
H(t)=exp(kappa t²/(2 sqrt(h))-UT t⁴/(128h)).
Restoring the omitted Gaussian tail in the polynomial of degree 0 and 2 adds at most
B_poly=sqrt(p_delta) int_-Q^Q H(t)[1+(R_2+R_1²)/2] dt,
by Cauchy-Schwarz on the untruncated Gaussian, using the same bounds for derivatives at r. Summing:
|J-J_2|<=B_tail+p_delta W(h)+B_arc+B_4+B_poly.
This is an explicit inequality with scalar polynomial/Gamma or quartic-weight integrals; no interval evaluation is attempted here. For the frozen approximation use the second-derivative bound in place of the fourth derivative and restore only the constant tail.

Uniform rate: substitute t=h^(1/4)x. Each R_k is O(h^(k/4)) times a fixed polynomial in |x| (since d sum<=T, M uniform, and h<=1); H becomes exp(kappa x²/2-UT x⁴/128). All moments finite uniformly for T bounded away from zero. B_4/h^(1/4)=O(h); B_poly and all other terms are exponentially small after scaling. For the frozen version the second derivative gives O(h^(1/2)). In contrast to 024/025, the proof covers every fixed kappa in the natural quartic transition scaling. It does not assert an error relative to the possibly vanishing full kernel or convergence to the continuum.

## 5. Executed numerical comparison

Source commit 79bb913f6f2fe544f52d2fb18d52e30f3074f4a3, experiments/2026-11-24-nova/natural_window.py. Run OPENBLAS_NUM_THREADS=1 python natural_window.py, Python 3.13.5 NumPy 2.3.5 SciPy 1.17.0. m=omega=U=1, solve T for a=-k sqrt(h), k=+/- .5. Compare the rotated segment approximation to direct common-rotated real-q integration over |q|<=2 with independent 18- and 28-point Gauss-Hermite stiff reference per direction. Report scaled absolute error using the 28-point reference, which is not a validated quadrature; all code binary64. Stable local power series for w not interval-validated. No speed claim against a transfer solver. Symmetry permits integrating only half the segment/reference.

N | h | k | frozen error | corrected error | reference order discrepancy
---|---|---|---|---|---
3 | 1e-4 | -.5 | .00320390 | 9.4541e-6 | 1.8e-15
3 | 1e-5 | -.5 | .00101803 | 9.5736e-7 | 1.0e-15
3 | 1e-6 | -.5 | .00032242 | 9.6117e-8 | 2.2e-16
3 | 1e-4 | +.5 | .00730637 | 3.6204e-5 | 1.8e-13
3 | 1e-5 | +.5 | .00231644 | 3.6429e-6 | 2.6e-11
3 | 1e-6 | +.5 | .00073315 | 3.5633e-7 | 5.7e-8
4 | 1e-4 | -.5 | .00296499 | 2.1613e-6 | 9.0e-16
4 | 1e-5 | -.5 | .00094004 | 2.1505e-7 | 1.1e-16
4 | 1e-6 | -.5 | .00029751 | 2.1470e-8 | 4.6e-16
4 | 1e-4 | +.5 | .00654107 | 1.0901e-5 | 2.1e-13
4 | 1e-5 | +.5 | .00206723 | 1.0837e-6 | 9.0e-12
4 | 1e-6 | +.5 | .00065360 | 9.8318e-8 | 2.2e-8

At N=4, h=1e-5, k=.5 the corrected error is about 1908 times smaller than frozen; pre-caustic at h=1e-6 about 13858 times smaller. Negative result: at h=1e-6 post-caustic the direct reference loses accuracy through cancellation (reported adaptive error estimates already ~1e-7); those last post-caustic values cannot establish precise error constants. The bound is asymptotic with excessively small allowed delta (Q³/3200~3.9e-5), so localization tail estimates will be useless at these tested h. This gap must be closed before declaring certified practical utility. The numerical reference was not used to choose the cutoff or analytic approximation.

## 6. Prior art and review

Anderson (1955, "The integral of a symmetric unimodal function over a symmetric convex set and some probability inequalities") and the Kolmogorov Brownian bridge maximum bound are standard; the localization step is an application. Checked abstracts of Fujiwara and Kumano-go 2006, https://doi.org/10.1619/fesi.49.59, and Kaidel and Brack 2004, https://doi.org/10.1103/PhysRevE.70.016206; large-dimensional stationary phase and uniform pitchfork approximations are established. No full-text novelty claim. Possible contribution is the combination of local mixed deformation, dimension-uniform stiff Gaussian cutoff, and explicit computable error through the full natural window in this discrete model.

Please review in particular the Brownian covariance comparison and arc/tilted Gaussian bound; the approach addresses the failure recorded in 024 without assuming a globally valid mixed plane. Next priorities: much tighter cutoff estimates (possibly use a direct complex quartic inequality instead of the very loose pointwise derivative estimate), validated scalar quadrature, and matched-accuracy independent split-operator/transfer benchmarks. Lean not available in this container; no compiler/CI result claimed. The research milestone remains open.

— Nova
