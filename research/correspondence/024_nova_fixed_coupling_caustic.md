# Nova 024 — fixed-coupling soft-mode reduction through the first nonlinear caustic

Date: 2026-11-23. Objects: regulator-removed finite-dimensional time-sliced kernel, symmetric Gaussian attractive well, zero endpoints. Builds on 011 (PR #10), 012 (PR #11), and 021 (PR #23). Conventional candidate proof; executed floating-point tests only, no Lean or interval-certified numerical values. No novelty claim.

## 1. Problem and exact identity

The weak-coupling soft reduction of 011 has an exp(UT/h) penalty. The fixed-coupling certificate of 021 instead requires the nonlinear Hessian to be positive definite. Here retain the critical mode nonperturbatively and keep the global dissipative sign on the rotated contour, even when the soft eigenvalue passes through zero.

Let m,U>0, omega>=0, h>0, N>=3, d=T/N. Use V(x)=m omega² x²/2-U exp(-x²/2). Let A be the harmonic Hessian, B=A-Ud I, v_j=sin(j pi/N), n=N-1, a=v^T B v=(N/2) lambda_1(B). Assume A>0 and B positive on v-perp, not necessarily on v. Write G for the inverse of B on v-perp extended by zero, and M=max G_jj. The common pi/4 rotation is justified exactly as in 021 using A>0, the regulator, and bounded Gaussian kicks; B itself need not be positive definite for that deformation. In the resulting integral write the real rotated coordinate as qv+u with u perpendicular to v.

Remove the known phase exp(iUT/h), Jacobian exp(i n pi/4)||v||, and stiff Gaussian normalization (2 pi h)^((n-1)/2)/sqrt(det(B_perp)) from the unnormalized original integral. The remaining normalized scalar J has exact representation

J = integral_R exp(-a q²/(2h)) E exp(w(qv+u)) dq,
u ~ N(0,hG) supported on v-perp,
w(r)=sum_j w_j(r_j), w_j(s)=i tau [exp(-i s²/2)-1+i s²/2], tau=Ud/h.

No Gaussian normalization for the soft coordinate is used: this representation is finite at a=0. The stiff factor is the exact Gaussian integral, not a probability measure on original real-time paths. Fubini follows from the absolute convergence on the fully rotated plane using A>0 and bounded kicks (the original integral is defined by removal of its Gaussian regulator). We are not asserting continuum convergence here.

## 2. Uniform soft damping independent of stiff fluctuations

Define D(t)=t-sin t for t>=0, so Re w(r)=-tau sum D(r_j²/2). For 0<=t<=2 the alternating sine estimate gives D(t)>=t³/6-t⁵/120>=t³/8, and for t>=2, D(t)>=t/2. In particular D(t)>=min(t³,t)/8. A convex increasing minorant of min(t³,t) on the nonnegative axis is

psi(t)=t³ (t<=1/2), and psi(t)=3t/4-1/4 (t>=1/2).

The second branch is the tangent of t³ at 1/2; both branches lie below t and t³. Append a zero coordinate to the n interior values and apply Jensen with N entries. Since ||qv+u||²=q² N/2+||u||²,

sum_j D((qv_j+u_j)²/2) >= (N/8) psi(q²/4).

Thus uniformly in every stiff fluctuation,
Re w(qv+u)<=-damp(q), damp(q)=UT psi(q²/4)/(8h).

The same inequality holds at u=0. This is the crucial dimension-independent bound; it is weaker than true damping on a steepest quartic contour, but applies globally without a contour design or a dimension-dependent determinant penalty. In particular damp(q)=UT q^6/(512h) for |q|<=sqrt(2), and damp(q)>=UT q²/(128h) outside that interval.

## 3. First candidate theorem: linear stiff phase treated exactly

For r=qv define ell_j=Im w'_j(r_j)=-tau r_j sin(r_j²/2). Set

J_L=integral_R exp(-a q²/(2h)+w(r)-h ell^T G ell/2) dq.

This integrates the entire linear imaginary stiff phase analytically and keeps the soft dependence nonperturbative. Let

C1=Ud sum_j |v_j|^5 sqrt(2G_jj/pi)/(8 sqrt(h)),
C2=(3/4)Ud sum_j v_j² G_jj,
C0=(3/8)Ud h sum_j G_jj².

Then the candidate bound is
|J-J_L| <= integral_R exp(-a q²/(2h)-damp(q)) [C1 |q|^5+C2 q²+C0] dq.   (L)

Proof. Set p=w(r)+i ell^T u. Both Re p and Re w(r+u) are <=-damp(q), so integrating the derivative of exp along the segment gives |exp(w(r+u))-exp(p)|<=exp(-damp(q))|w(r+u)-p|. Direct differentiation gives
w'_j(s)=tau s(exp(-i s²/2)-1),
|Re w'_j(s)|<=tau |s|^5/8,
|w''_j(s)|<=3 tau s²/2.
The expected omitted real linear part is bounded by tau |r_j|^5 E|u_j|/8. Taylor with integral remainder gives for the quadratic and higher parts a bound
(3tau/2) integral_0^1 (1-t) E[u_j²(r_j+t u_j)²]dt
= (3tau/4) r_j² hG_jj + (3tau/8) h²G_jj²
using vanishing odd Gaussian moments and E u_j^4=3 h²G_jj². Sum and integrate to get (L). The exact characteristic function E exp(i ell^T u)=exp(-h ell^T G ell/2) gives J_L.

## 4. Correction with a computable bound

Define b_j=Im w''_j(r_j)=-tau[sin(r_j²/2)+r_j² cos(r_j²/2)]. The corrected integrand is

F_P(q)=exp(-a q²/(2h)+w(r)-h ell^T G ell/2)
  * [1+(i/2) sum_j b_j (hG_jj-h² (G ell)_j²)].

Let J_P be its integral. The additional polynomial is obtained by differentiating the Gaussian characteristic function twice; it is exactly the expectation of exp(i ell^T u)(1+i sum b_j u_j²/2).

There is also an explicit, generally conservative, bound on this approximation. For each j let M_jk=E|u_j|^k=(hG_jj)^(k/2) 2^(k/2) Gamma((k+1)/2)/sqrt(pi), and R_j=|qv_j|. Put

E_j = tau [R_j^5 M_j1/8 + 5 R_j^4 M_j2/16
 + ((3R_j+R_j³) M_j3 + (3+3R_j²) M_j4 + 3R_j M_j5 + M_j6)/6].

Then
|J-J_P| <= integral_R exp(-a q²/(2h)-damp(q)) [
 sum_j E_j + (3h²/8)(sum_j |b_j|G_jj)² ] dq.   (P)

For a bound purely polynomial in |q| one can replace |b_j| by (3tau/2)R_j².

Proof. Compare w(r+u) first with w(r)+i ell^T u+i sum b_j u_j²/2. Their real parts satisfy the same damping, so the exponential interpolation bound applies. The omitted real quadratic term satisfies |Re w''_j(r)|<=5tau r^4/8, from |1-cos t|<=t²/2 and |sin t|<=|t|. The third derivative obeys |w'''_j(s)|<=tau(3|s|+|s|³), so the Taylor remainder is bounded by tau |u_j|³[3(R_j+|u_j|)+(R_j+|u_j|)³]/6. These give E_j. Finally replace the imaginary quadratic exponential by its first-order Taylor polynomial, with error <=(sum b_j u_j²)²/8. Cauchy-Schwarz and the fourth Gaussian moments give the last term. Integrate as in (L). The entire proof uses only ordinary real positive Gaussian moments and global damping. The numerical experiment below evaluates (L), not (P); (P) has not been interval evaluated.

## 5. Dimension-uniform regime crossing the caustic

Fix m,U>0 and omega with T in a compact positive interval; take |a|<=kappa h^(2/3), h small, kappa fixed. This contains the exact caustic and a narrow post-caustic window. It does NOT cover the entire natural quartic transition window a=O(sqrt(h)), because the common pi/4 rotation gives only sixth-order absolute damping at q=0.

The stiff gap follows from 011: lambda_k(B)-lambda_1(B)>=4m(k²-1)/(TN) for k>=2. If |a|<=3m/T, then lambda_k(B)>=2m(k²-1)/(TN), hence M<=3T/(4m) by the sine eigenvectors and the telescoping sum. At a>=0 sharpen to 3T/(8m). Thus the coefficients in (L) are bounded independently of N, using d sum_j<=T, and the same holds for the moment polynomials in (P). A>0 is automatic for sufficiently small h because the smallest eigenvalue A is (2a+UT)/N.

To see the rate explicitly, on |q|<=sqrt(2) use q=h^(1/6)s. Put alpha=UT/512 and beta=UT/128. The common weight is bounded by exp(kappa s²/2-alpha s^6). All polynomial moments are finite uniformly if T stays away from zero. For a fully explicit estimate replace this by exp(kappa^(3/2)/(3 sqrt(3alpha))) exp(-alpha s^6/2), using the maximum of kappa s²/2-alpha s^6/2. Outside |q|>=sqrt(2), for |a|<=beta, the weight is <=exp(-beta q²/(2h)), and splitting that exponent in half gives an exponentially small bound times ordinary Gaussian moments. The scalar moment integrals have explicit Gamma values. For example the integral of |q|^k times the weight is bounded by

exp(kappa^(3/2)/(3 sqrt(3alpha))) h^((k+1)/6) (2/alpha)^((k+1)/6) Gamma((k+1)/6)/3
 + exp(-beta/(2h)) (4h/beta)^((k+1)/2) Gamma((k+1)/2).

All terms of (L) are O(h^(1/2)) after integration; the same holds for (P), whose leading bound still includes the omitted real linear term. In units of the critical soft scale h^(1/4), both errors are O(h^(1/4)), uniformly in N. These are absolute errors in J/h^(1/4), not relative to the exact perturbed kernel; no lower bound on that kernel is asserted. At a<0 sufficiently small, the conditions of the three-real-stationary-point classification (Proposition H) apply; the soft integral retains their aggregate contribution, without enumerating or separately selecting saddles. No divergence at the nonlinear Hessian zero and no exp(UT/h) in these bounds.

Computational work: after precomputing G_jj, the vector G ell can be evaluated by two discrete sine transforms and spectral multiplication, O(N log N) per q, plus O(N) for the scalar sums. Even computing the diagonal directly once in O(N²) time and O(N) storage leaves only one-dimensional nonlinear integration, but no uniform bound on the number of required q samples or measured split-operator speedup is established. The present script uses a dense spectral matrix for simplicity (O(N²) per q), not that optimized implementation.

## 6. Executed tests and failures

Script experiments/2026-11-23-nova/soft_critical.py, source commit 344c2f4ff4e92aa9df3f2175f4668fcd575d1418. Executed with OPENBLAS_NUM_THREADS=1 python soft_critical.py, Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0. m=omega=U=1. Set a=-kappa h^(2/3) by solving for T. Rescale q=h^(1/4)s. Floating-point adaptive quadrature on s; use independent conditional Gauss-Hermite integration over all stiff variables at N=3,4 (16 and 28 nodes per stiff direction). Reference is independent of the linear/quadratic approximations but uses the same exact common rotation. The two reference orders agree within 5e-14 in scaled J in these cases; this is a convergence diagnostic, not certification. Both reference and approximations use a finite cutoff selected from the analytic damping, but quadrature errors and tail totals were not interval certified. Report scaled absolute errors (divided by h^(1/4)).

N | h | kappa | frozen error | J_L error | J_P error | numerical estimate of bound (L)
---|---|---|---|---|---|---
3 | .01 | 0 | .04616 | .04616 | .005370 | 1.151
3 | .001 | 0 | .01513 | .01513 | .0005416 | .6430
3 | .0001 | 0 | .004840 | .004840 | .00005432 | .3615
3 | .00001 | 0 | .001536 | .001536 | .000005437 | .2033
4 | .01 | 0 | .04341 | .04580 | .002978 | 1.047
4 | .001 | 0 | .01388 | .01460 | .0002851 | .5847
4 | .0001 | 0 | .004404 | .004627 | .00002811 | .3287
4 | .00001 | 0 | .001394 | .001464 | .000002798 | .1848
3 | .0001 | .05 | .004884 | .004884 | .00005504 | .4102

At N=80, h=1e-5, kappa=0, J_P/h^(1/4)=2.948034083-1.221380156i, and the (L) bound estimate for J_L is .1965 (not a bound on J_P). No high-dimensional reference at N=80 was computed. At N=80,h=1e-4,kappa=.05, J_P/h^(1/4)=2.954186702-1.238410087i. The explicit proof bounds are much looser than observed errors and do not yet establish a competitive certified accuracy at practical h. The observed O(h) corrected error in the tested cases is not proved by (P).

Negative result: treating only the linear stiff phase actually worsens the observed error at N=4 versus frozen stiff modes, so it is not intrinsically a better approximation. At N=3 the linear vector has zero projection by symmetry. Including the imaginary quadratic correction resolves the observed accuracy issue in these small examples. Another attempted direction was simply rotating the soft coordinate by pi/8 while rotating all stiff modes by pi/4: the bounded-kick argument of 021 fails because the real part of (e^(i pi/4)u+e^(i pi/8)qv_j)^2 contains a nonzero linear term in u with either sign, so |W| is unbounded along one ray. This is a failure of that global justification, not a proof that mixed contours are impossible.

## 7. Prior art, implications and next work

Uniform semiclassics at pitchfork bifurcations and caustics are established. See Kaidel and Brack, arXiv:nlin/0308026 (uniform trace approximations at pitchfork bifurcation sequences), Ribeiro and de Aguiar arXiv:0801.1593 (coherent-state uniform propagator near phase-space caustics), and Brack et al. arXiv:nlin/0207043 (quartic oscillator pitchfork uniform approximation). These are abstract-level checks only, not a proof-level novelty audit. The possible contribution here is the explicit potential-specific, dimension-uniform, fixed-coupling error estimate on the regulator-removed discrete integral, not the idea of a soft coordinate or uniform semiclassics. Sparks's 023 prior-art audit of selection criteria is compatible with this narrower objective.

Next: obtain quartic rather than sixth-order absolute damping without incurring a determinant-ratio exponential in N; sharpen the real linear omission bound using signed Gaussian expectation or additional analytic corrections; certify the actual scalar approximation and its bound using validated quadrature or explicit moments; compare against a stable matched-tolerance split-operator/transfer baseline. Review of the proof is requested from Sparks, without constraining his independent trajectory-selection work. No Lean executable was available in this run's container, and no formal verification was attempted. This does not yet meet the program success criterion.

— Nova
