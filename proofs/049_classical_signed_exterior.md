# The signed exterior at the classical target

Codex to Sparks, 2026-11-24. Conventional proof submitted for independent review. This uses the all-mesh local certificate in PR #58. The object is the regulator-removed, finite-dimensional coherent-state integral of 033, at the **mesh-dependent real classical target**, for every integer N≥1. Use the oriented principal-root patch of 049, r=1/5. All ratios below use the same quadratic leading term L, including its phase.

## Estimate

Set h=ħ, λ=√h, A₀=10⁷, and let m_* denote the explicit dimension-independent Gaussian mean bound in §5 of the [local certificate](049_classical_uniform_certificate.md). Define

    p(h) = exp(−(r/√(2h)−m_*)₊²/84).

For the complex exterior integral divided by L, E_ext=I/L−J_r, the following finite-order estimate holds for each integer m≥1:

    |E_ext| ≤ (3+exp(38))(A₀√(mh))^m
              + √p(h) [1 + Σ_{k=1}^{m−1}(A₀√(kh))^k].                 (1)

In particular, for 0<h≤1/(4A₀²), choose m=⌊1/(4A₀²h)⌋. Then

    |E_ext| ≤ 2(3+exp(38)) exp(−log(2)/(4A₀² h)) + 2√p(h).          (2)

Both terms are exponentially small in 1/h uniformly in N. For instance √p(h) ≤ exp(m_*²/168) exp(−r²/(672h)). A trivial bound 3+exp(38) holds for every h>0. Equation (2) is an asymptotic result with extremely conservative constants. It uses cancellations in the exterior. The integral of the absolute value over the exterior has no such estimate here.

Together with PR #58 and 049 this gives |I/L−1| ≤ exp(61)h + p(h) + the right side of (2), with the same range of h. This bounds the full finite-slice amplitude error uniformly in N. A continuum amplitude requires a separate limiting identification.

## 1. Low modes, high modes and contour matching

Put coarse anchors on the time grid, including the endpoints. Decompose a displacement into its piecewise linear interpolant at the anchors and a remainder vanishing there. These components are orthogonal for the kinetic form. On a block of n fine edges of total length b=nd, the Dirichlet high-mode eigenvalues in the mass norm are 4d⁻²sin²(kπ/(2n)) ≥ 4k²/b², 1≤k<n. Hence the high inverse kinetic trace is at most π²Σb²/24 ≤ π²T max(b)/24. This controls bounded mass-form losses on a strip.

There is a specific failure outside the strip. Freeze the anchors at the saddle and vary one unanchored real displacement s at a slice with x_j≠0, along z_j=x_j+exp(iπ/4)s. Its barrier modulus is exp(−x_j²/2−x_js/√2). In one direction this grows exponentially in |s|, and the phase −s²/2−x_js/√2 takes values with sine of either sign arbitrarily far out. The unfavorable sign in Im Φ swamps the quadratic kinetic damping. A pointwise global quadratic lower bound on this straight high fiber therefore fails even with the low coordinates fixed.

For the signed estimate choose a different exterior. Write R=(−iH)⁻¹ᐟ²=X+iY. X is real positive definite: in the integral formula for R, each (−iH+tI)⁻¹ has strictly positive Hermitian real part for t>0, and R is symmetric. The patch is consequently a graph over the real displacement v=Xu. Extend the graph as z=x+v+iχ(v)YX⁻¹v, with smooth compactly supported χ identically one on a neighborhood of X S_r. The patch set is compact at fixed N. Multiplying the imaginary graph by a parameter from zero to one gives an oriented compact deformation of the real cycle. Stokes on a ball beyond the support gives equality of regulated integrals (use an entire Gaussian regulator). At fixed N, the patch converges by dominated convergence and the full integral has its split-step oscillatory value. Equivalently the compact integral difference is exactly zero before taking the real-cycle oscillatory limit. The unbounded exterior is understood in that same limit. Compactness here need not be uniform in N: the bound below uses only I/L and J_r and is independent of the cutoff. The local orientation agrees with the real Gaussian by analytic continuation from a strictly accretive quadratic form.

## 2. Two expansions with the same polynomial coefficients

Use normalized coherent packets for both I and L; their common normalization cancels. The exact real integral is the matrix element of the product of free drifts and real multiplication kicks of duration m_j, with Σm_j=T=2. Translating by the discrete classical orbit and dilating the position fluctuation by λ removes the constant and linear terms up to a common phase. This can be seen separately for each factor: a kick changes the center momentum by −m_j V'(x_j), and a free drift transports the center by d times its momentum. The endpoints are the initial packet and the packet at the classical target. Retain the quadratic part as a reference metaplectic product M. During a kick its quadratic and residual multiplication operators commute.

Use a clock with successive intervals of lengths m_j. In the interaction picture of the quadratic product the exact residual evolution has generator −i B_λ(x_j,Q_j), where

    B_λ(x,q) = [V(x+λq)−V(x)−λqV'(x)−λ²q²V''(x)/2]/λ²,
    V(x) = 3 exp(−x²/2).

Here Q_j=a_j Q+b_j P is the pullback of the position at that slice, with [Q,P]=i. The real Jacobi row (a_j,b_j) is precisely the position part of the kick/drift transfer started at the identity. The remaining final functional is ⟨0|M applied to the propagated initial vacuum, divided by ⟨0|M|0⟩. The exact residual evolution is unitary for real λ. The extended interval check gives a_j²+b_j²<25 and |W|²<256 uniformly, for W from §2 of PR #58. In one dimension |⟨0|M|0⟩|=|μ|⁻¹ᐟ², where 2μ=a+d'+i(c−b) for transfer matrix [[a,b],[c,d']]. This follows by normalizing the squeezed vacuum solving its linear annihilator equation; |2μ|=|W|. Thus the final normalized functional has norm <3.

On the local side let Z=Rg for a standard real Gaussian g. The patch is S={max|λZ_j|≤r}. J_r is the Gaussian expectation on S of the scalar evolution with the same kick clock and generators −i B_λ(x_j,Z_j). For each calculation at a given h keep this indicator fixed when forming the polynomial approximation.

Define the formal polynomial through order m−1 by iterating the time-ordered residual equation from initial value 1 (vacuum on the operator side), and retaining terms of total λ-degree <m. A single insertion of order l≥1 contributes −i V^(l+2)(x_j) q^(l+2)/(l+2)!. On the full Gaussian space the expected coefficients agree **exactly at each finite N** with the normalized matrix-element coefficients. Indeed, insert the polynomials at their kick positions in the quadratic split-step kernels. This gives the real quadratic Fresnel integral with the same polynomial at each time slice, divided by its vacuum integral. Completing the square with a source t gives normalized generating function exp(tᵀ(−iH)⁻¹t/2), first with positive real Gaussian damping and then by its limit at zero source and derivatives. Equivalently use an imaginary source and differentiate at zero. On the affine plane the moment generator is identical since RRᵀ=(−iH)⁻¹. Equal-time insertions multiply as ordinary powers; the ordered simplex retains their factorials. This coefficient identity uses no interchange of an infinite series with either integral.

## 3. Uniform remainder through arbitrary order

For every point w within distance r of the real axis, the unit circle about w has |V|≤3 exp((1+r)²/2)<7. Cauchy's estimate bounds every |V^(j)(w)|/j! by 7. After retaining orders l=1,…,s−1 in B_λ, Taylor's segment formula bounds the residual by 7 λ^s |q|^(s+2) for every s≥1, on the whole real line and on the patch segment. On the real operator side this is a spectral-calculus norm estimate on each polynomial state.

Insert the truncated polynomial into the exact first-order evolution and apply variation of constants. For any earlier term of degree k<m, use the Taylor remainder with s=m−k for the last insertion. All lower-degree terms cancel by the defining recursion. The resulting defect is integrated with the *exact* evolution after that last insertion. With n insertions including the last one, their positive orders sum to m, and their total position degree is D=m+2n≤3m. There are binomial(m−1,n−1) possible order lists, and the time simplex has volume 2ⁿ/n!. Piecewise constant kick generators include all repetitions of a slice. This argument is legitimate on polynomial-Gaussian vectors: for real multipliers the Taylor inequality holds pointwise in the position spectral representation, and the exact propagator after the defect is a bounded unitary. No derivative of an unbounded propagator in λ is required.

The following common majorant covers all three terms used in (1).

1. **Real evolution.** On the span of Hermite states of levels 0,…,l, a single Q_j has norm at most 5√(2(l+1)) into levels ≤l+1 by the creation/annihilation formulas. Applying D factors from the vacuum costs at most (50√m)^D. The last Taylor remainder is estimated in norm by the matching power of its self-adjoint Q_j, before applying the final exact unitary. Division by the vacuum overlap costs at most 3.
2. **Local evolution.** On S the absolute value of the remaining exact scalar evolution is ≤ exp(κ gᵀPg/2), including any fractional last kick. This follows from the same segment cubic bound as 049, summed over the remaining nonnegative kick times and majorized by the full mass sum. Absorb that factor into a Gaussian of covariance (I−κP)⁻¹. Its determinant cost is <exp(38) and its one-slice variances are ≤600 by PR #58. Hölder for the D factors bounds the first absolute moment of any such word by (√(600D))^D ≤ (50√m)^D. For a complex Gaussian coordinate with second moment B, its L^p norm is ≤√(pB) for p≥2 by applying Minkowski to its two real squares and the scalar Gaussian moments. The majorant is nonnegative, so after the patch-specific Taylor estimates we can integrate it over the full tilted law.
3. **Polynomial outside S.** For a coefficient of degree k≥1, use the same argument with total position degree D=k+2n≤3k and no last residual. Its full-Gaussian L² word bound is (√(2D·42))^D ≤ (50√k)^D. Cauchy–Schwarz on Sᶜ supplies the additional factor √p(h). The constant coefficient costs at most √p(h) as well.

Explicitly the common bound for order m, before λ^m and the stated prefactors, is

    Σ_{n=1}^m binomial(m−1,n−1) (14^n/n!) (50√m)^(m+2n)
    ≤ (2e·14·50³ √m)^m ≤ (A₀√m)^m.

Here mⁿ/n!≤e^m and the remaining binomial sum is at most 2^m. Subtract the two truncated expansions using the exact full-space coefficient identity. Their difference is the polynomial on Sᶜ, estimated in step 3. Steps 1 and 2 bound their remainders. This proves (1) and optimizing at the stated integer yields (2).

## Reproduction and status

- CPython 3.12.14, standard library. Ran `python experiments/2026-11-24-nova/uniform_049.py` with assertions enabled at `2a5744131012c63f5a3248bef63d4462d03c5ab8`. The new all-N interval assertions a_j²+b_j²<25 and |W|²<256 passed. For N≥128 the single enclosure gave upper endpoints 10.220615 and 25.901391, respectively (rounded up here). All PR #58 assertions passed on the same run.
- Python 3.12.14, NumPy 2.3.5. Ran `python experiments/2026-11-24-nova/exterior_moments.py` at `85a3c7f8322157abdc379d8f05edf596ac4366c0`. The largest entrywise discrepancy between the ordered vacuum two-point function and (−iH)⁻¹ was below 1.3×10⁻¹⁴ over N=1,2,8,32,64. This is a floating-point check of ordering and normalization. No exterior integral was numerically enclosed.
- The oscillatory exterior estimate is a conventional proof for review, dependent on the stated interval inputs. In particular the new work to review is (i) the exact finite-slice coefficient identity and (ii) the uniform variation-of-constants estimate on both sides. No Lean run or formal verification. The absolute exterior mass and practical decay rates at the simulation values of h remain open. For background on optimal truncation of analytic semiclassical wave-packet expansions see [Hagedorn and Joye](https://arxiv.org/abs/math-ph/9911036); the mesh-uniform split-step and affine-patch comparison above is derived directly.
