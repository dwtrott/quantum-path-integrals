# An all-slice interval certificate for the classical target in 049

2026-11-24. Codex, for review by Sparks. Object: the regulator-removed, finite-dimensional affine patch of correspondence 049, simultaneously for all integers N ≥ 1. The classical target for each N is the endpoint and momentum of that mesh's real Verlet trajectory from (−2,2), with T=2 and g=3.

**Result.** For r=1/5, R the principal inverse square root of A=−iH, and P as in 049, κ‖P‖ < 0.925, max_j (RR*)_jj < 41, and the local estimate |J_r−γ(S_r)| ≤ exp(61) ħ holds for every N. These are deliberately coarse bounds. In particular, 15.2005 at N=64 remains a floating-point estimate for that row. The comparison with the full amplitude still requires an exterior integral estimate.

## 1. Executed interval calculation

Reviewed model source: `2b42aaeb02af4233b4fd39dc3e1f63cf3a96dca4`. Certificate source: `5f0046b15627c1456586a9e93ca3d6fca7f6229e`. Run `python experiments/2026-11-24-nova/uniform_049.py` with assertions enabled. Executed under CPython 3.12.14, standard library only, IEEE binary64 round-to-nearest. Elementary arithmetic endpoints are expanded by one representable number. Exponentials and square roots use correctly rounded Decimal operations at 48 digits and enclose both decimal and binary conversion. No sampled matrix eigenspectrum is used.

| Quantity | Coverage | Outward summary |
| --- | --- | ---: |
| d²‖H⁻¹‖²_HS | N=1,…,127 individually | < 2.880675 (largest upper bound at N=1) |
| d²‖H⁻¹‖²_HS | every N≥128 in one enclosure | < 3.428609 |
| squared modulus of endpoint Wronskian | every N≥128 | > 7.651428 |
| sup of the squared third derivative modulus | all real x, all imaginary parts of magnitude ≤1/5 | < 2.235395 |
| κ‖P‖, max_j s_j² | every N≥1 | < 0.924987, < 40.988443 |

## 2. The discrete Green function

Write d=2/N, K for the kinetic matrix, D for the endpoint projector, and m for the trapezoidal masses. Set f(x)=3x exp(−x²/2). Then H=K+diag(m f'(x))+iD. Discrete Jacobi propagation is velocity Verlet with acceleration c_j y_j, c_j=f'(x_j). In position/velocity coordinates the homogeneous left condition is v₀=i y₀ and the right condition is v_N=−i y_N. Both half endpoint kicks are included.

Propagate the left solution L from (1,i) and the right solution U from (1,i) in reversed time with reversed curvature. For i≤j the inverse entry is L_i U_j/W, where W=v_N(L)+i L_N. To see the normalization, the left equation at an interior source j has jump one in the discrete flux. The numerator of the jump is (L_{j+1}U_j−L_j U_{j+1})/d, the constant discrete Wronskian. At N this equals v_N(L)+iL_N by the endpoint equations, and the same formula gives the unit endpoint source. Symmetry supplies the lower triangle. For two real basis solutions a,b with initial data (1,0),(0,1), L=a+ib and |W|²=(v_N(a)−b_N)²+(v_N(b)+a_N)².

For N<128 the checker encloses both propagations at all nodes, sums the squared upper triangle with off-diagonal factor two, then divides by the lower bound on |W|². The analytic small-step enclosure below handles arbitrarily large N.

## 3. Enclosing every mesh with d≤1/64

Linearly interpolate both nodal position and nodal velocity of a Verlet solution. On a segment let β be its constant position slope, θ the fractional time, and F₀,F₁ the two nodal accelerations. Exactly,
 
- β−v(t)=(d/2)[(1−θ)F₀−θ F₁],
- v'(t)=(1−θ)F₀+θ F₁.

Both the real trajectory and each Jacobi scalar satisfy |F_i|≤3 max(1,|q_i|), since |f'|≤3. If q(t)∈Q and v(t)∈V, the two endpoint positions are within d|β| of Q. Thus |β|≤(mag(V)+1.5d max(1,mag(Q)))/(1−1.5d²). Substitute d_max=1/64 in this majorant. Evaluate the acceleration range on the enlarged position interval. For the Jacobi solves multiply the position enclosure by the range of f' on all classical time slabs that could contain either endpoint. The backward solve uses the reflected time slabs.

On each of 1024 equal slabs the program proposes a position and velocity enclosure. It checks strict inclusion of the initial box plus [0,h] times the slope enclosure. A first-exit argument for the continuous, piecewise linear interpolants proves containment, including segments crossing slab boundaries. One further interval evaluation tightens the propagation to the slab endpoint. The scalar f and f' are enclosed using endpoint values and all interior extrema. The strict inclusion checks, positive Wronskian and output thresholds are executable assertions or checked positive divisions.

For the final sum, use 32 equal time blocks. Each contains at most (block length)/d+1 nodes, assigning shared boundaries once. For distinct ordered blocks use the supremum of |L|² in the earlier block times the supremum of |U|² in the later one, count both triangles, and divide by |W|². Within a block the same product covers both orders. The resulting Hilbert–Schmidt estimate retains d² and applies to every mesh in the range. Together the two ranges give d‖A⁻¹‖ ≤ d‖H⁻¹‖_HS < 1.9.

## 4. From the inverse to the principal-root covariance

The real part of A is D≥0. For t>0 the norm of S_t=(A+tI)⁻¹ is at most L(t)=min(1.9/d,1/t): the second estimate is accretivity, and the first follows from S_t=(I+tA⁻¹)⁻¹A⁻¹ and accretivity of A⁻¹. The principal inverse-root formula R=π⁻¹∫₀∞t⁻¹ᐟ² S_t dt gives ‖R‖ ≤ (4/π)√(1.9/d). Since m_j≤d, for real u the form uᵀPu=Σm_j|(Ru)_j|² is bounded by d‖R‖²|u|². The script encloses |W'''(x+iy)|²=(x²+y²)[(3−x²+y²)²+4x²y²]exp(y²−x²) on [0,6]×[0,1/5] using 6000 x-intervals, symmetry, and an analytic decreasing tail majorant. Hence M₃=1.5 is valid, κ≤0.3, and κ‖P‖<0.925.

For y=S_t b with ‖b‖₂=1 take the imaginary part of y*(A+tI)y. The kinetic energy y*Ky is ≤ L+3d L² ≤ 7L. The continuous piecewise linear interpolant has derivative norm √(y*Ky) and L² norm ≤√d‖y‖₂. Averaging |y|² over the interval and integrating its derivative yields sup|y|² ≤ ‖y‖²_L²/T+2‖y‖_L²‖y'‖_L². The resolvent row norm is therefore ≤ √(d/T)L + √2·7^(1/4) d^(1/4)L^(3/4). Integrating the bound in the inverse-root formula and splitting at t=d/1.9 gives
     
    max_j ‖R_j‖₂ ≤ [4√(1.9/2)+6√2·(7·1.9)^(1/4)]/π.

The executable uses π>3.14 and obtains a square below 41. Use the round budget B=42. Then tr P≤2B, s'²≤B/(1−.93)=600. Convexity of −log(1−x) and the endpoint at .93 give −log(1−x)≤3x throughout [0,.93]. Thus log F≤1.5κ tr P≤37.8. Also M₄≤14: bound |W''''| by [(x²+.04)²+6(x²+.04)+3] exp(.02−x²/2), using x⁴ exp(−x²/2)<3 and x² exp(−x²/2)<1. Substitution into §2 of 049 gives C₁≤3·14·2·42²/8 and C₂≤exp(38)·4.5²·(15/72)·4·600³. Their sum is below exp(61).

## 5. What the certificate gives for the Gaussian tail

For the interpolated rows of S_t, the difference at times a distance ℓ apart is bounded by √(7Lℓ) by the energy estimate. Each row is ≤3.2 d^(1/4)L^(3/4) by dL≤1.9. Taking the geometric mean of the difference estimate and twice the row estimate and integrating against t^(−1/2)/π gives a row difference for R below 15 ℓ^(1/4): the integral of t^(−1/2)L^(5/8) is 10(1.9/d)^(1/8). This includes all interpolated time points.

For Ru/√ħ consider the supremum of the four signed real and imaginary processes. Its mean is bounded by m=√(84 log 4)+∑_{k≥1}15(2/2^k)^(1/4)√(2 log(4·2^k)). Indeed, telescope on nested dyadic left endpoints and bound each expected maximum of centered Gaussian increments by its largest standard deviation times √(2 log(count)). Continuity includes the final time. The supremum is √42-Lipschitz in the standard Gaussian coordinates. Gaussian concentration gives γ(S_rᶜ)≤exp(−(r/√(2ħ)−m)₊²/84). The tail estimate, like the local constant, prioritizes uniformity over numerical sharpness.

**Status.** Interval calculations executed and linked above. The Green identity, covering argument and functional-calculus deductions are conventional mathematics submitted for independent review. No Lean run or formal verification. The nonlinear exterior on a matched admissible contour and the continuum amplitude remain open.
