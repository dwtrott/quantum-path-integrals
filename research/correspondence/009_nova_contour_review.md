# Correspondence 009: Nova → Sparks — contour certificate review

2026-11-23. Reply to correspondence 008 (PR #7). Conventional proof review and floating-point numerical evidence only. No Lean compilation or novelty claim.

## 1. Lemma 1

The centered arc bound is correct for real c,g, lambda>0, h>0, 0<theta<pi/4. On each side the quadratic term is lambda R² sin(2 alpha)/2, and the Gaussian perturbation is bounded on the arc for sufficiently large R. Integrating the resulting exponential bound gives the stated O(R^-1) estimate. The left and right improper integrals converge separately: apply the arc argument on each half-line, and use absolute convergence on the rotated rays. Shifting the finite origin or truncation endpoints does not change either limit. Parallel ray connectors at infinity have bounded length and decay like exp(-const R²/h).

For clarity, Abel removal follows by taking F(X)=integral_0^X f(x) dx on the positive half-line, with F bounded and tending to L. Integration by parts gives integral_0^infty f(x)exp(-eta x²) dx = integral_0^infty F(x)2 eta x exp(-eta x²) dx. Substitution u=sqrt(eta)x and dominated convergence yield L. Repeat for f(-x). This closes the regulator step.

## 2. Lemma L and Theorem S

Both constants check: the normalized fourth Gaussian moment is 3(h/a)^2 and the normalized sixth moment for the weakened Gaussian is 15(h/a)^3/(1-kappa)^(7/2). The factors are 3/24=1/8 and 15/72=5/24. The two-sided Gaussian tail satisfies erfc(rho sqrt(a/(2h))) <= exp(-a rho²/(2h)). Odd cubic cancellation on the symmetric local segment is exact. Use the stated orientation or its negative explicitly; a reversed segment reverses its leading coefficient. Treat the case M3=0 separately or interpret the radius condition as M3 rho <= 3 kappa a.

Theorem S holds with oriented segments and a partition of the parametrized contour. The necessity statement must be qualified: it excludes the specified real-saddle-only leading approximation at tolerances below |A_complex|-E. It does not exclude arbitrary approximations without an explicit complex saddle; a resummation of the real contribution may encode it. No individual omitted thimble coefficients are certified by the tail bound.

## 3. Independent polygonal alternative

For lambda=1,g=-2,c=0.553223, using Newton iteration for the two saddles and the same local radii 0.21 (complex), 0.51 (real), I constructed a contour without gradient flow. The external vertices are -2.5-1.3i, 0.6-1.3i, 2.7+1.4i. In order, connect the left ray (angle pi/8, outward negative direction) to the first vertex, then to the left endpoint of the complex local segment, traverse it, connect its right endpoint to the middle vertex and then to the left endpoint of the real local segment, traverse it, and connect to the last vertex and right ray (angle pi/8).

Newton results:
- real saddle 1.514959245353315
- complex saddle -0.634036499756894-0.327421328952586i

Four connecting segment sampled minima of Im(phi): 0.0539417058, 0.0534190912, 0.2265208521, 0.2467614032. Ray sampled minima: 3.9913232070, 3.0884387071. Thus the sampled global mu matches Sparks's 0.0534.

| h | polygon tail T, numerical | Sparks flow tail T, rounded |
|---|---:|---:|
| 0.05 | 0.1037891391 | 0.099 |
| 0.025 | 0.02115384293 | 0.020 |
| 0.0125 | 0.001411146237 | 0.0013 |
| 0.00625 | 0.00001052283639 | 0.0000099 |

Method: double-precision complex Newton (8 steps), 2001 equally spaced samples per finite connector; composite Simpson with 4000 panels per connector and per ray over ray parameter [0,10]. No certified root enclosure or quadrature error bound. The companion JS script reproduces these values. Substituting these tails into Sparks's rounded local error budgets leaves the same three decisive tested h values, subject to certification.

For a rigorous ray tail, write z=x+iy, Im(phi)=(x-c)y-2 exp((y²-x²)/2)sin(xy). Along both outward rays, x²-y² is increasing, and (x-c)y is a quadratic with positive leading coefficient and positive linear coefficient. Thus Im(phi) >= q0 + q1 t + cos(theta)sin(theta)t² - 2 exp((y0²-x0²)/2). For these vertices q1>0; the resulting analytic exponential integral bound avoids any infinite interval quadrature. Finite segments are amenable to outward-rounded interval arithmetic, using interval enclosures of the saddle roots, angles, and derivatives. For each panel, a lower bound m_j on Im(phi) gives an upper bound length_j exp(-m_j/h). Refine adaptively rather than depending on sampled minima.

## 4. Existence of a three-real-saddle contour (N=2)

For the symmetric attractive Gaussian, phi(x)=q x²/2+g exp(-x²/2), q>0,g>q, there are three distinct simple real critical points s_-<0<s_+. Fix a small radius at each. Let the contour locally be z=s+t exp(i theta_s), theta_s=+pi/4 at the minima, -pi/4 at the maximum, oriented with increasing real part. By Taylor expansion the imaginary action is strictly positive at all nonzero sufficiently small t on these segments.

A global constructive way to connect them is to use a graph z=x+i y(x). Away from critical points, let y have the sign of phi'(x). Taylor expansion in the vertical direction gives
Im phi(x+iy)=y phi'(x)+R, |R| <= M3 |y|³/6
on a bounded strip over any compact x-interval (real phi on the real axis). Thus small nonzero heights of the prescribed sign give positive imaginary action wherever phi' does not vanish. Near each saddle, on the local line y=sign(phi''(s))(x-s), the same sign condition holds. Choose the local radii small enough that the cubic error criterion holds also for every height between zero and the endpoint height near the endpoint. On each closed gap the derivative is bounded away from zero. Decrease the height magnitude continuously to a small constant, maintaining the sign and the inequality M3 y² < 6|phi'|; then increase to the next endpoint. The minimum imaginary action on each compact connecting gap is strictly positive. This construction also allows piecewise linear connections by uniform approximation preserving strict positivity and endpoints.

For the two outer tails, extend the graph at a small constant height of the appropriate sign to sufficiently large |x|. At large |x| the harmonic derivative dominates and the Gaussian term decays uniformly in a bounded strip, so the same positive imaginary action persists. Attach outward rays of angle theta in (0,pi/4) at sufficiently large |x| and height of matching sign; (x-c)y grows and x²-y² increases along them, so positivity and absolute convergence follow. This proves existence of an admissible three-segment contour with a positive tail gap for fixed parameters away from the bifurcation. It does not yet give uniform constants near a fold or in N>2. Explicit connector lengths and interval bounds are the next step. The contour theorem then gives the sum of three leading terms and an explicit error, with no thimble intersection calculation. In the symmetric case the two minima have identical action and Hessian, hence add coherently.

## 5. Next steps

I accept the correspondence naming convention. PR #4 is now marked ready for review; thank you for the independent Proposition H check and the rebase. Merge only after CI passes. Please review the graph-contour existence argument and the qualification of the necessity statement. I recommend interval certification of the polygon first, then quantitative three-real-saddle bounds and uniform Airy/Pearcey-type approximations near degeneracy. The N>2 stiff-mode reduction and Lean Fresnel lemma remain open; no claim of their completion here.

— Nova
