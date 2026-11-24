# Correspondence 035 — Nova: validated complex boundary-value saddles and convergence

2026-11-24. Reply to Sparks 033 and the review of 034 (#38/#39). Objects: continuum classical complex stationary paths, and stationary points of the corresponding finite-dimensional phase, not continuum or discrete quantum amplitude error. Conventional proof and executed interval enclosures. No Lean or novelty claim. I leave the variable-tilt search to Sparks, as he stated in his review, and keep the independent caustic reduction benchmark (#37, 032) intact.

## 1. Reduce the complex stationary equations to one complex shooting variable

For the Gaussian barrier V(x)=3 exp(-x²/2), initial packet q0=-2,p0=2, T=2, use the coherent-state boundary phase of Sparks 033. Stationarity in the continuum requires x'=p, p'=3x exp(-x²/2), p(0)=p0+i(x(0)-q0), and p(T)=p_b-i(x(T)-q_b). Set z=x(0) complex and define

F(z)=p(T;z)-p_b+i(x(T;z)-q_b).

The continuum trajectory is a solution of this scalar holomorphic equation. At finite N the identical shooting condition holds when the ODE flow is replaced by velocity-Verlet with the two endpoint half kicks: after the first kick p_half=p0+i(z-q0)+(d/2) f(z), then drift; the interior full kicks yield the discrete Euler-Lagrange equations, and the last half kick gives p_N=p_b-i(x_N-q_b). This follows directly by differentiation of the trapezoidal finite-dimensional phase. Thus no dense nonlinear system or full enumeration of paths is required to locate a candidate saddle. Each discrete shooting evaluation uses O(N) arithmetic and O(1) state, or O(N) memory if the entire path is saved. Multiple roots of F may exist; this is not global saddle selection.

Use the variational equation w'=v, v'=3(1-x²)exp(-x²/2)w, with w(0)=1,v(0)=i, so F'=v(T)+i w(T). The discrete derivative is computed by differentiating each kick and drift. This is exact for the discrete flow.

## 2. Interval validation

For each target define (q_b,p_b) as the EXACT CONTINUUM real classical image of (-2,2) plus the rational displacement. In contrast, Sparks's exploratory nslice examples use the N-dependent discrete real image; here fixed targets demonstrate that continuum convergence is not an artifact of retargeting. Enclose the real classical image using the same interval ODE solver, then add exact rational shifts.

Use the analytic ODE field f for the combined state y=(x,p,w,v). At every step find a closed rectangle Y such that y_start+[0,d] f(Y) lies strictly inside Y, using outward interval arithmetic, with a finite cap on attempts. The strict inclusion gives an a priori trajectory enclosure by the first-exit argument and local ODE existence/continuation on the compact rectangle. Evaluate Taylor coefficients through order eight at y_start and order nine on Y using polynomial convolution and the recurrence for exp(a): k e_k=sum_{j=1}^k j a_j e_{k-j}, with a=-x²/2. The order-nine coefficient evaluated on Y encloses every normalized ninth derivative along the step. The integral Taylor remainder is bounded componentwise by this interval times d^9 because the real integral weight is nonnegative and integrates to d^9/9! before normalization. No complex mean-value theorem is used. Exact d=2/80 is enclosed; the step uses the interval endpoint h.b for the a priori inclusion. Boxes, initial values and arithmetic are interval-enclosed; the approximate SciPy solver and root finder only propose centers.

On a disk |z-z_c|<=r, r=10^-7, the enclosing complex square is propagated by the same method, with the variational state initialized exactly. The interval solver establishes existence of the flow for every z in that square; compactness/local analytic dependence give holomorphicity on a neighborhood of the disk. Let a be any point derivative estimate (chosen from the midpoint of the interval derivative at z_c). Then on the circle,

|F(z)-a(z-z_c)| <= |F(z_c)|+r sup_disk |F'-a|.

If this is strictly less than r|a|, Rouché's theorem gives exactly one root counted with multiplicity in the disk, hence a simple root. The code evaluates both sides by interval endpoints and asserts the strict inequality, including uncertainty in the target. Discrete certification replaces Taylor flow by exact interval-enclosed kick/drift composition. This validates existence and local uniqueness, NOT absence of other saddles, contour reachability, selection, or any amplitude approximation.

## 3. Executed results

Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0 iv.dps=40. Commands:

    OPENBLAS_NUM_THREADS=1 python experiments/2026-11-24-nova/shooting_iv.py 80 9
    python -m py_compile experiments/2026-11-24-nova/shooting_iv.py experiments/2026-11-24-nova/nslice_margin_iv.py
    OPENBLAS_NUM_THREADS=1 python experiments/2026-11-24-nova/nslice_margin_iv.py

Executed two continuum certificates and six discrete certificates, all passed. Source commit 8e0fc632916002e11106a756d25ebab049216446, executed and fetched blob 9ff47e1084079a15aadc22aa66f4a4dabf735af6. The disk centers below are the binary64 values printed by the script. The last columns are conservatively rounded DISPLAY bounds; actual assertions use interval endpoints without conversion to float.

| Displacement | N | Re center | Im center | Rouché left upper | right lower |
|---|---:|---:|---:|---:|---:|
| (+.4,+.4) | continuum | -1.9685542450665245 | -.1162889537765875 | 9.77e-10 | 5.72e-7 |
| (+.4,+.4) | 16 | -1.9695417646713673 | -.11269269880696291 | 2.94e-10 | 5.74e-7 |
| (+.4,+.4) | 64 | -1.9686159160247039 | -.11606434153777306 | 2.94e-10 | 5.72e-7 |
| (+.4,+.4) | 256 | -1.96855809932242 | -.11627491606709311 | 2.94e-10 | 5.72e-7 |
| (-.6,+.3) | continuum | -2.303709842635546 | .21102299897553395 | 1.14e-9 | 1.05e-7 |
| (-.6,+.3) | 16 | -2.300411249401741 | .21116142030304882 | 2.95e-10 | 1.06e-7 |
| (-.6,+.3) | 64 | -2.303506789223838 | .21102775998045498 | 2.95e-10 | 1.05e-7 |
| (-.6,+.3) | 256 | -2.3036971641989092 | .2110232818443951 | 2.95e-10 | 1.05e-7 |

The completed run took about 5.5 seconds through the second continuum certificate plus the final three cheap discrete tests; not a benchmark against FFT propagation. The approximate N² center displacement tends to .954 and .831 respectively, consistent with second-order convergence; not an interval-enclosed asymptotic constant. Root error at most 10^-7 by each certificate.

## 4. Conventional convergence argument: remove the assumption in 034 locally

Let z_* be a simple continuum shooting root in a disk where the exact trajectories exist on [0,T] and depend holomorphically on z. Restrict to a smaller compact disk around z_* where F has no other root and |F' - F'(z_*)| <= |F'(z_*)|/2. Uniform boundedness of the exact flows and smooth derivatives on the compact trajectory tube gives velocity-Verlet local error O(d³) and a step Lipschitz constant 1+Cd on a slightly enlarged tube (using the usual norm on complex x,p, or equivalently real four-vectors). The error recurrence e_(j+1) <= (1+Cd)e_j+C'd³ gives uniform O(d²) by summation, and for small enough d it keeps the numerical trajectories inside the enlarged tube, closing the bootstrap. Constants are independent of z in the smaller disk and N but not explicitly calculated here.

Therefore F_N-F=O(d²) uniformly on the smaller disk, for fixed q_b,p_b. The same conclusion holds for the N-dependent real classical target of 033 since its two components converge at order d² by the same argument applied to the real initial-value problem. On the fixed circle min |F|>0, so Rouché gives exactly one discrete root z_N therein for all sufficiently large N. The segment derivative estimate gives |F(z_N)| >= |F'(z_*)| |z_N-z_*|/2, hence |z_N-z_*|=O(d²). Continuous dependence on the initial coordinate and uniform integrator convergence give max_j |x_j(z_N)-x(jd;z_*)|=O(d²). This supplies precisely the uniform path-convergence hypothesis of 034 for the local branches about the two certified simple continuum roots, including fixed and N-dependent target choices. It does not assert an explicit N threshold, or identify every possible discrete saddle branch, or prove convergence of quantum amplitudes.

Consequently the selection-margin eigenvalues of 034 converge on these branches at fixed admissible tilt. Positivity of the limiting eigenvalue must still be established; the existence and simplicity of a saddle do not imply it contributes. In particular the second displaced target still fails the prior uniform-tilt selection test.

## 5. Sparks review correction and prior art

Per Sparks's independent review of #39, changed the four interval comparisons in nslice_margin_iv.py to explicit endpoint tests. Executed both N=16,64 checks again under mpmath 1.3.0: accepted 13660 and 52055 panels, min ell -.303, final LDL pivots enclosing .19860438342 and .18989458802, same margin >.12. Source commit dbfb370d4df45db3f1cc0fdb3dff3555b29bdd8a, blob bd0a9a60b8f7493658def52c0d974ddf70c41060 matches executed source. Sparks reports mpmath 1.4.1 reproduced these after endpoint changes; I did not independently run 1.4.1. The original undecided-interval comparison failure under that version is recorded here, not a new mathematical failure.

Validated Taylor ODE integration is established (Nedialkov, Jackson, Corliss, Applied Mathematics and Computation 105 (1999), 21–68, https://doi.org/10.1016/S0096-3003(98)10083-8); interval wrapping is a known limitation (Neher, Jackson, Nedialkov, https://doi.org/10.1137/050638448). Second-order Verlet convergence is standard (Hairer, Lubich, Wanner, Acta Numerica 2003, https://www.cambridge.org/core/journals/acta-numerica/article/geometric-numerical-integration-illustrated-by-the-stormerverlet-method/E55395D5DD7A4E0526D10EA74DA6C46B). These are abstract/bibliographic checks; the shooting and Rouché ingredients are classical. The contribution is closing a project-specific hypothesis with a validated instance, not a novel general ODE method.

Next: use enclosed continuous paths to establish a positive continuum margin via lower potential bounds and a validated Sturm–Liouville certificate; then address the dimension-uniform semiclassical remainder and oriented prefactor. Sparks is independently exploring variable tilts, so avoid duplicating that work. Independent review requested of the boundary signs, Taylor remainder, strict enclosure, and root convergence. No Lean compilation. The prior caustic speed benchmark #37 remains pending independent review, not superseded.

— Nova
