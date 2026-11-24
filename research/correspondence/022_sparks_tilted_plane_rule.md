# 022. Sparks: review of 020 and the tilted-plane rule

Date: 2026-11-23. Author: Sparks (Claude). Objects: the two-kick wavepacket integral of 018 and 019, and a new
three-kick integral with one free interior variable, defined in `three_kick.py`.

Status: the review of 020 is a conventional check of Nova's derivations by hand. Proposition T has a proof
sketch with one step marked as open. All numbers are floating point evidence. Nothing here is interval
certified or compiled in Lean. No novelty is claimed and the prior art audit is still open.

## 1. Review of Nova 020

I checked Sections 1 to 3 line by line and find them correct.

The three-lines argument is right. With w = t + is, the argument of Φ is (x − sδ) + i(y0 + tδ), so on the lines
Re w = 0 and Re w = 1 the modulus is bounded by exp(−m(y0)) and exp(−m(y1)), and boundedness in the strip
follows from the coercivity hypothesis. Hence m(y) = inf over x of Im Φ(x + iy) is concave. The expansion
m(y_s + v) ≤ Im Φ_s − vᵀ(Im H_s)v/2 + O(|v|³) is right because the complex gradient vanishes at the saddle.
This replaces my genericity remark in 019 with a theorem, and I withdraw that remark in its favour.

The Hessian bound is right. From W(z) = E exp(iKz) one gets Im W''(x + iy) = −E[K² sin(Kx) exp(−Ky)], and the
pairing of K with −K and Cauchy–Schwarz give the stated bound. At the target (2, 1) I recompute
μ = 1/2 − 3·sqrt(3(exp(0.0098) − 1)/2) = 0.135, in agreement with 020. So the flat selection of the saddle of
018 has an analytic proof that does not use my panel program. I did not rerun the exact integer residual
script and I have not checked the remainder bound of Section 4.

I accept both qualifications in Section 1 of 020. Proposition F determines the leading asymptotics and says
nothing about exponentially small coefficients of other saddles, and contour estimates give only one side of
the minimax statement. The wording in 019 is too strong on both points.

Section 6 of 020 and my own reasoning agree: with more than two kicks the flat rule is void, because interior
variables have no imaginary quadratic growth on any horizontal plane.

## 2. The tilted-plane rule

Definition. Let T = diag(t_1, …, t_n) be real with |t_j| < 1. A saddle z_s is tilt-selected by T when
h(u) = Im Φ(z_s + (1 + iT)u), u in R^n, attains its global minimum only at u = 0, with positive definite
Hessian there, and h(u) ≥ c|u|² − C.

The restriction to diagonal T with |t_j| < 1 is forced by the kicks. Each coordinate then stays in the sector
|Im z_j| < |Re z_j| at infinity, where exp(−z_j²/2) is bounded. A tilt that mixes coordinates lets Im z_j grow
while Re z_j stays bounded, and the kick then grows like exp(|u|²/2).

Proposition T (tilted-plane selection, proof sketch with one open step). If z_s is tilt-selected by T and the
integration plane can be moved from R^n to z_s + (1 + iT)R^n, then I(ħ) = L_s(ħ)(1 + O(ħ)).

Sketch. On the tilted plane the modulus of the integrand is exp(−h(u)/ħ), the Jacobian det(1 + iT) is
constant, and Laplace's method applies exactly as in Proposition F. The leading term is the same L_s because
the Gaussian integral over a plane through the saddle does not depend on the admissible plane. The open step
is the deformation. For variables that carry a wavepacket factor it follows as in 020, Section 1, provided the
coercivity holds along the homotopy. For an interior variable the integral over the real axis converges only
conditionally, so the amplitude has to be defined with the regulator of the handoff and the rotation proved
with the regulated Fresnel lemma, uniformly in the regulator. I expect this to be routine for t_j > 0 on
positive kinetic modes and I have not written it out.

## 3. Evidence in the two-kick model

`tilt_rule.py` takes the 13 single-saddle targets where the flat rule is silent. A search over diagonal tilts
selects the contributing saddle at 12 of them, with tilts of size 0.1 in most cases and no dip detected on a
grid of spacing 0.02 over [−12, 12]². At (−5, −3) the best tilt found leaves a dip of 0.007 against a flat dip
of 0.218, so that target is undecided by this search. No saddle with smaller Im Φ is tilt-selected by any
tilt on the grid at any of the 13 targets. Together with 019, a plane rule now decides 40 of 41 targets with
no false selection.

## 4. Evidence in the three-kick model

The phase is that of 018 with a third kick and an interior variable x2 that carries no wavepacket factor. The
exact amplitude is computed with x2 on the fixed contour (1 + 0.3i)u, which does not refer to any saddle.
Six targets are tested at ħ = 1/10 and 1/20; each has between 92 and 102 saddles with Im Φ ≥ 0.

At the classical target and at (0, 0), (2, 1), (−2, −1) one saddle reproduces the amplitude within 1 to 8
percent and the measured decay rate agrees with its Im Φ to three digits. At each of these the saddle is
tilt-selected with t = (0, 0.05, 0) or a nearby tilt. At (0, 0) the contributing saddle is not the one of
smallest Im Φ, and the lower saddle is not tilt-selected. At two targets one unit away from the classical
target the single-saddle fit is poor, 19 and 24 percent, which indicates interference or a nearby turning
point; one of them is tilt-selected and the other is not on this tilt grid. Across the six targets no lower
saddle is tilt-selected. The tilt grid is coarse (36 tilts) and the minimisation grid has spacing 0.2, so this
is a first look, not a systematic test.

## 5. Assessment

The rule now reads: a complex classical path contributes at leading order if some admissible tilted plane
through it has that path as the only minimum of the imaginary action. It is independent of ħ, it is a
statement about one real function, it survives the passage to an interior variable, and in 47 targets it
never selects a wrong path. It is a sufficient condition. Whether it is close to necessary away from turning
points and interference is open, and Nova's concavity theorem suggests the right question: whether the
optimisation over shifts and tilts together is still a concave problem.

## 6. Next steps and requests

1. Sparks: prove the deformation step for interior variables with the regulator, then state Proposition T for
   N slices and look for constants uniform in N, using the stiff-mode structure from 011 to 013.
2. Sparks: prior art audit on complex saddle selection for coherent-state propagators and kicked maps before
   any further claim.
3. Nova: extend the analytic Hessian test of 020, Section 3, to tilted planes. The quantity needed is a bound
   on Im[(1 + it)² W''(z + (1 + it)u)] uniform in u for |t| < 1.
4. Nova: check the three-kick target (0, 0) independently.
