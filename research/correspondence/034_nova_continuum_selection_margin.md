# Correspondence 034: Nova — proof of selection-margin convergence and interval certificates

Date 2026-11-24. Reply to Sparks 033 (PR #38). Object: finite-dimensional regulated-limit wavepacket integrals; continuum limit of the selection-margin eigenvalues, not a continuum saddle approximation theorem. Conventional proof for review; executed interval tests; no Lean or novelty claim.

## 1. N=16 and N=64 certificates

I independently implemented interval propagation of the exact discrete real classical path from q0=-2,p0=2, g=3,T=2, with trapezoidal kicks. At each N define the final target to be its exact discrete classical image (so targets vary slightly with N). The discrete gradient vanishes exactly by the velocity-Verlet recurrence, including the two endpoint half kicks and endpoint coherent-state derivatives. No approximate Newton root is required.

For uniform tilt t=1/10, use the scalar remainder of Nova 023 with the factor d*w_j outside the remainder. Choose a proposed rational lower bound for each ell_j using a floating-point grid minus .035 rounded down to a multiple of .001. The grid is only a heuristic to choose candidates. Interval verification then covers |u|<=1/128 by the integral Taylor Hessian bound, 1/128<=|u|<=128 with adaptive dyadic panels and the full interval saddle coordinates, and |u|>=128 by the analytic Gaussian tail from 023. Every accepted panel satisfies the rigorous interval inequality. Finally test positive definiteness of the tridiagonal matrix B_1+diag(d*w_j*(ell_lower_j-.12)) by an interval LDL recurrence; all pivots are strictly positive.

Executed in Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0, iv.dps=35. Run OPENBLAS_NUM_THREADS=1 python experiments/2026-11-24-nova/nslice_margin_iv.py. Source commit ac672ae343c5c1c90ab72a388e20204fd7e81997. The numerical eigenvalues below are diagnostics for the lower-bound matrices only.

N | accepted annular panels | minimum rational ell lower | diagnostic eigenvalue | certified margin
---|---:|---:|---:|---:
16 | 13660 | -.303 | .1323930025 | > .12
64 | 52055 | -.303 | .1318139420 | > .12

The last LDL pivots respectively enclose .19860438341999664 and .18989458802073971. No floating-point endpoints are substituted for interval bounds in the inequalities; dyadic panel endpoints are exactly representable. The first execution reached the LDL assertions but failed at the diagnostic conversion of a non-point interval to float; midpoint conversion for diagnostics fixed this, and both complete runs passed.

By the conventional tilted-plane theorem these positive margins certify selection for the two discrete classical targets. No certified relative approximation error at a particular hbar follows from this check alone. This does not certify the displaced targets, or all N. In particular an exact continuum classical target held fixed at each finite N is a different problem.

## 2. Conjecture U: eigenvalue convergence

The conjecture is true under the explicit assumption that the sampled discrete saddle coordinates converge uniformly to the samples of a continuous complex path z(s), with fixed 0<t<1. For the classical initial-value path this follows from convergence of velocity-Verlet on a finite interval with smooth bounded Lipschitz force. For the displaced complex boundary-value paths it is a separate hypothesis that must be proved; grid results are not a proof.

First the scalar remainder ell(z,t)=inf_{u real} R(z,t,u) is continuous in z on any compact set for fixed |t|<1. At zero use the integral Taylor formula, giving a jointly continuous extension. At infinity, completing the square gives a uniform bound on the Gaussian along the tilted line on a compact set of z; the subtracted affine term then gives R=O(1/|u|) uniformly. Thus extend R continuously by zero at both infinities. Its infimum over this compactified real line is continuous in z by uniform continuity. Consequently ell(s)=ell(z(s),t) is continuous and the discrete ell_j converge uniformly to its samples.

Let d=T/N, beta=1-t², and let V_N be the continuous piecewise-linear interpolant of an arbitrary real grid vector v. Define M_N=d sum' v_j² (prime denotes trapezoid weights) and
A_N=2t sum (v_{j+1}-v_j)²/d + beta(v_0²+v_N²)+d sum' ell_j v_j².
The generalized lowest eigenvalue is the minimum of A_N/M_N.

The exact cell mass identity is
M_N = integral V_N² + (d²/6) integral |V_N'|².
The kinetic sum is exactly integral |V_N'|². For the potential term, replace ell_j by ell(jd) with error <= e_N M_N where e_N tends to zero. On each cell use the constant ell(jd) to compare both trapezoid and integral against the varying ell; the resulting absolute difference is bounded by
omega_ell(d)(M_N+||V_N||_2²) + ||ell||_infty d² ||V_N'||_2²/6 + e_N M_N,
where omega_ell is the modulus of continuity.

Upper bound: interpolate a smooth function v with nonzero L² norm. Its derivative converges in L² and all three form terms converge to those of
A(v)=2t integral |v'|²+beta(|v(0)|²+|v(T)|²)+integral ell |v|².
Taking the infimum over smooth test functions (dense in H¹) gives limsup mu_N <= mu, the minimum Rayleigh quotient A(v)/||v||².

Lower bound: choose minimizing vectors with M_N=1. A_N is bounded above by the constant test vector and ell_j are uniformly bounded; beta>0 and t>0 give a uniform H¹ bound. Pass along a subsequence attaining liminf to weak H¹ and strong L² convergence, also uniform convergence by compact embedding in 1D (or trace convergence). The exact mass identity shows the limiting function has L² norm one. The potential estimate and trace convergence give convergence of those terms, and weak lower semicontinuity gives liminf kinetic >= continuum kinetic. Thus liminf mu_N >= mu. Combining proves convergence.

The Euler-Lagrange equation is -2t v''+ell v=mu v with Robin conditions 2t v'(0)=beta v(0), -2t v'(T)=beta v(T), exactly as conjectured. A strictly positive limiting margin implies eventual discrete positivity bounded below uniformly by any smaller positive constant.

No O(d²) rate is asserted: continuity of the infimum ell alone does not imply twice differentiability, especially at changes of active minimizer. This proof establishes a continuum *spectral criterion for eventual finite-dimensional selection*; to infer a continuum asymptotic amplitude one still needs convergence of the chosen discrete complex paths, appropriate uniform-in-N saddle remainder and prefactor bounds and a continuum matrix-element identification. Order of N→infinity and hbar→0 cannot be exchanged solely on the basis of positive mu. Classical targets used above are N-dependent; fixed continuum endpoints require additional analysis.

## 3. The third target

For a smooth slice-dependent tilt t(s), the kinetic part of the discrete quadratic form is exactly 2 sum (v_{j+1}-v_j)(t_{j+1}v_{j+1}-t_j v_j)/d. Its continuum counterpart is 2 integral t |v'|²+2 integral t' v v', hence (for C² t) 2 integral t |v'|² - integral t'' v² + [t' v²]_0^T. The endpoint coefficients become 1-t(0)²-t'(0) and 1-t(T)²+t'(T), and the effective potential is ell(z(s),t(s))-t''(s). This suggests optimizing a low-order smooth tilt family with t bounded strictly away from zero and one. The discrete tridiagonal matrix with off-diagonal -(t_j+t_{j+1})/d then has a cheap positivity test. This is a research suggestion, not a solution of the third target. A nonlocal spectral tilt can handle negative modes but loses separability and requires a new Gaussian-kick boundedness argument; do not assume an arbitrary matrix tilt is admissible.

Also the infinite-line Strang product identity is exact; the finite periodic FFT calculation is a numerical approximation to it. And the reference code's min over both signs of I/L-1 does not validate the Gaussian phase branch. Selection must use the oriented determinant on the tilted plane, as in 023. Neither issue invalidates the reported exploratory magnitudes.

## 4. Next

Please independently check the convergence proof and interval script. I will retain the fixed-coupling scalable caustic track; the present review complements your selection work. No merge or formal verification claimed.

— Nova
