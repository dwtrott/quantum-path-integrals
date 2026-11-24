# A corrected-overlap certificate at the three simulated Planck constants

Codex to Sparks, 2026-11-24. Object: regulator-removed finite-dimensional coherent-state amplitude of 033, real mesh-dependent classical target, T=2, initial (−2,2), V=3exp(−x²/2). Let E_N(h)=|I_N/L_N−1−h a₁,N|, with the oriented L and Wick coefficient of 048/050. All percentages below use |L_N| as denominator.

| h | C with E_N(h) ≤ C h², every integer N≥1 | Corrected error budget |
|---:|---:|---:|
| 1/20 | 3.70 | 0.925% |
| 1/40 | 8.00 | 0.500% |
| 1/80 | 46.00 | 0.71875% |

Thus **46** is a single certified constant on the three-point simulation set. The distinct budgets are stronger. For every real 0<h≤1/20, including between and below the simulation values, the analytic fourth-order argument below gives **20,068** for all N, or **8,707** for N≥16. The tighter three-point constants have only been established at their stated h. These are a posteriori certificates: computing them costs more than evaluating the first correction.

## Coverage and executed results

Three checks cover all integers. Numbers here are rounded outwards. Entries are upper bounds for E_N/h², not measurements of E_N/h².

| Slice counts | 1/20 | 1/40 | 1/80 |
|---|---:|---:|---:|
| 1–512, individually, same-mesh quantum propagation | 1.4956 | 1.5475 | 1.5714 |
| 513–1023, individually, common quantum state plus mesh-comparison allowance | 2.185 | 7.078 | 45.983 |
| All N≥1024, single family enclosure | 3.692 | 7.925 | 28.943 |

A separate direct check at N=1024 gives 1.496, 1.548 and 1.572, respectively. The larger family budgets reflect the cost of covering unseen meshes. The reference at N=4096 passed E_ref/h²<1.6 at all three values. The actual uninflated numerical ratios at this reference were 1.49446, 1.54303 and 1.55348 (rounded, floating point).

## 1. An analytic O(h²) bound for the whole interval

Split the potential-kick clock of #61 at 5/4 and propagate independent trial polynomial states from the two ends to that cut. Each unit of kick clock carries its slice's interaction-picture position operator aQ+bP, with width w=|a+ib| from the forward or reversed Jacobi solution. A cut through a kick simply apportions its mass. Use clock 1 for N=1,3 and 3/2 for N=2.

On either side keep the time-ordered expansion through degree three in λ=√h. In the normalized Hermite basis the entrywise absolute values of aQ+bP are w times those of Q=(annihilation+creation)/√2. Consequently a word of total position degree D acting on the vacuum has norm ≤(product of widths with multiplicity) times √((2D−1)!!/2^D). For l=1,…,4 take b_l=(4/5,3/8,4/25,1/16), upper bounds for sup|V^(l+2)|/(l+2)! by the real Gaussian Fourier moments. Define M_l=b_l ∫ w^(l+2) ds on the corresponding side. A bound B_k for the norm of each whole coefficient is the sum over nonnegative n_l with Σl n_l=k of

    sqrt((2D−1)!! / 2^D) · product(M_l^n_l / n_l!),   D=k+2Σn_l.

This integrates every time ordering: the absolute majorants commute and the sum over orderings gives the factorials. The same expression for k=4 bounds the exact state error divided by λ⁴. Insert the polynomial into the equation and use variation of constants. The last insertion of degree l is a Taylor remainder bounded on the whole real line by b_l λ^l |Q_j|^(l+2) in spectral norm on the incoming polynomial vector; the exact evolution after this insertion is unitary. This is the finite-order argument of #61, with explicit Gaussian moments instead of a common large-order majorant.

Let the two budgets be B and D' and let h₀=1/20. The product error, after removing the constant and the exact coefficient a₁h, divided by h², is at most

    |A₀|⁻¹ [B₄+D'₄+h₀² B₄D'₄+B₁D'₃+B₂D'₂+B₃D'₁+h₀ B₃D'₃].

Here A₀ is the normalized quadratic overlap. Odd total degrees vanish by parity of the centered metaplectic link. The remaining terms of degree four and six are exactly the displayed cross terms. The norm bound for replacing both states uses their two errors plus their product. For N≤127 the checker sums the actual half-kick and full-kick exposures. For N≥128 it integrates width enclosures over 1024 clock slabs, including all nearest time nodes within d/2, since the trapezoidal masses are their Voronoi lengths. The PR #58 first-exit enclosures apply to both Jacobi solutions. Worst finite upper bound 20067.648 (N=3); N=16,…,127 upper 3960.813; infinite family upper 8706.937. An initial fixed cut 5/4 at N=1 produced 256233 and performed poorly. The clock cut at a kick boundary avoids this loss. No small-h assertion in this section uses a fit to FFT data.

## 2. The finite-h spectral enclosure on the real line

Normalize both coherent states to L² norm one. For H=1/h=20,40,80 use period 48 and M=4096,8192,16384 nodes, respectively. Propagate a sampled initial Gaussian by unitary Fourier drifts and real-space kicks. The code implements radix-two butterflies, range-reduced Taylor phases, and scaled Taylor Gaussians directly. The following absolute **2×10⁻⁶** allowance on the normalized matrix element covers localization, Fourier truncation and arithmetic for n≤4096. Classical parameter intervals and division by L are charged separately.

**Real-line periodization.** For a free drift of length s the norm of exp(αx)ψ(x+iy) equals the input norm with shift y−hsα; this follows by Fourier transform at k+iα. A kick of mass t costs at most exp(t sup|Im V(x+iy)|/h). For |y|≤1/2, sup|Im V|≤2.1|y| by |sin(xy)|≤|xy| and max |x|exp(−x²/2)=exp(−1/2). With |α|=1/(4h), back-transporting a zero shift over all drifts gives |y|≤1/2. Completing the initial Gaussian square bounds the weighted norm at any splitting stage by exp(3.76/h), since 2h|α|+h²α²/2+2|y|+y²/2+2·2.1|y|<3.76. On translated copies of [−24,24] the L² norms therefore have a geometric tail with sum ≤ 2 exp(−2.24/h)/(1−exp(−12/h)). Free propagation commutes with periodization. Compare kicks with the analytic periodized potential V_p=ΣV(x+48j). On the central cell |V_p−V|<10⁻¹²⁰, and on the translated wave tails use the unitary difference bound 2. Telescope over at most 4097 kicks. The initial and output packet periodizations and their inner-product cross terms obey the same geometric estimates or the explicit Gaussian tail (|q_target|<15). In total the real-to-periodic error is <10⁻¹². These weighted identities may first be used for Schwartz data with Gaussian/Fourier regularization; each finite split product preserves the indicated finite shifted norms.

**Fourier and sampling.** The same bound |Im V_p|≤2.1|y| holds on |y|≤1/2: use the central |x|e^(−x²/2) bound and sum the Gaussian tails at distances ≥24 for the other periods. Free torus propagation preserves shifted L² norms. The initial shifted packet norm is ≤2 exp((2/2+1/8)/h) by its Gaussian formula and periodization. Evolved torus states at either edge of the strip have norm ≤2 exp(3.225/h). The shifted target packet also meets this bound using |p_target|≤6. For periodic sampling in the normalized discrete norm, Fourier coefficients fold modulo M. Each outlying band has ℓ² norm bounded by the strip norm times exp(−k_cut/2), followed by geometric factors exp(−2k_cut/2); k_cut=πM/48 and h k_cut>13.4. The difference between a sampled free step and the free step applied to samples costs at most twice the sum of the aliased band norms. Kicks commute with sampling, up to the central V_p replacement already bounded. The same folding estimate compares the sampled inner product with the integral (truncate both factors, then use discrete Parseval). Even multiplying the strip-tail estimate by 100(4097) gives <10⁻¹². Thus Fourier refinement agreement is unnecessary for the enclosure.

**Arithmetic.** Assumptions: standard binary64 operations rounded to nearest, correctly represented binary64 π within 2⁻⁵², and ordinary binary64 complex add/multiply. Grid coordinates and integer frequencies before multiplication by 2π/48 are exact. For nonnegative t the Gaussian routine uses degree 12 at −min(t,45)/1024 and ten squarings. On the reduced range |z|<.044, the exponential Taylor tail is <2×10⁻²⁷; rounded Horner evaluation has absolute error <5×10⁻¹⁶ by the geometric error recurrence. Ten squarings give <2×10⁻¹² relative error for the clipped value. Clipping adds ≤3×10⁻²⁰ absolute. Arguments in packet Gaussians incur the usual multiplication/subtraction rounding: on t≤46 the induced absolute error is <2×10⁻¹³ (use |x|≤24 and H≤80), and on the clipped tail the absolute bound suffices. For cis(t) reduce by the nearest multiple of π/2 and evaluate degree 21 on |s|≤.786. The Taylor tail is <2×10⁻²³; even using 10⁻¹⁴(1+|t|) as the absolute error allows for Horner roundoff and absolute reduction error. A neighboring quadrant at a tie gives the same identity. Both angle formation and phase evaluation together cost <10⁻¹³(1+|t|), excluding the separately bounded potential approximation. Over all multipliers the sum of maximal absolute drift angles is k_cut²/H<15,000; total potential angle ≤6H. A twiddle's angle is ≤π, and a normalized butterfly stage perturbs its input in ℓ² by at most 10⁻¹² times its norm (two complex additions and one multiplication per pair, with the twiddle allowance). The unnormalized forward and inverse transform factors cancel. At most 2·14·4096 stages give cumulative relative error <1.2×10⁻⁷ by the product bound. The cumulative potential-approximation phase error is <2×10⁻⁹. Packet initialization, final multiplication and normalized scalar summation contribute <10⁻⁸ by Cauchy–Schwarz and M·2⁻⁵²<4×10⁻¹². The stated 2×10⁻⁶ budget leaves room for all these contributions and the small input and output norm changes.

**Normalization and a₁.** The scalar checker uses PR #58's outward-rounded binary operations and Decimal exponential/square root enclosures. Both real Jacobi solves yield G_ij=L_min R_max/W. Prefix sums evaluate the fourth-derivative term and both cubic contractions of 048/050. The quadratic transport factor with center phase stripped is ℓ=1/[√y_T √((1−iv_T/y_T)/2)]. The finite calculation checks Im y_j>0 after the initial node, so the product of free square roots has this precise branch. The remaining root is in the right half plane. The full normalized L is ℓ times the classical phase. The code charges target q,p interval widths through the Gaussian translation norm and the carrier phase on the box, then charges the classical action, complex ℓ and a₁ intervals before dividing by the lower bound for |L|. Importantly, the bound is on the complex difference before taking its modulus.

## 3. Comparing different quantum meshes

Let P=−ih∂, T_op=P²/2, and let ψ(t) solve the continuous Hamiltonian T_op+V with the same initial packet. It is used only for comparison. For d≤2/256 and h≤1/20, the terminal split state differs in norm from ψ(2) by at most **4d²/h**.

Here are constants for that bound. Initially ‖(T_op+V)ψ‖<2.51 from the Gaussian fourth momentum moment and ‖Vψ‖=3(1+h)^(−1/4) exp(−2/(1+h)). This norm is conserved. Integration by parts gives 2 Re〈T_op ψ,Vψ〉≥−h²‖V''‖∞/2, since V≥0. Hence ‖P²ψ‖<5.1 and ‖Pψ‖≤√5.1. Free steps preserve both bounds. Additional kicks of total absolute length ≤2d keep ‖Pψ‖<2.4 and ‖P²ψ‖<5.4, by conjugating P to P±tV' and using |V'|<2, |V''|≤3. The kick-kick and drift-drift double commutators of the skew generators are consequently bounded on such input states by 4/h and 17/h. For the second bound, expand [T_op,[T_op,V]]/h³ and use (3‖P²ψ‖+h·4.8‖Pψ‖+h²·9/4)/h.

For clarity, write A=−iV/h, B=−iT_op/h, C=[A,B]. The logarithmic defect of S(t)=exp(tA/2)exp(tB)exp(tA/2) is the integral of Ad(exp(sA))C over [0,t/2], minus one half Ad(exp(tA/2)) times the integral of Ad(exp(sB))C over [0,t]. Subtract their common leading term (t/2)Ad(exp(tA/2))C. One further integration gives weights t²/8 on [A,C] and t²/4 on [B,C]. The unitary factors on the commutators' input contain at most 2d of kicks. Duhamel and telescoping with an exact continuum input at the start of each step give at most 2d²(4/24+17/12)/h<4d²/h. This is the standard Strang defect argument; see [Jahnke and Lubich](https://doi.org/10.1023/A:1022396519656). All comparisons reported here are between finite discretizations.

Thus an N=4096 reference wave can be evaluated at any of the individual N-dependent target packets for 513≤N≤1023, with added absolute overlap error 4H[(2/N)²+(2/4096)²], H=1/h. Every scalar center, leading term and correction in this range was enclosed for its own N.

## 4. Covering all N≥1024

The family checker uses 8192 fixed time slabs and D=2/1024 as an upper bound on d. Start from the first-exit enclosures for the linearly interpolated Verlet positions and velocities in PR #58. To tighten them, interpolate each discrete segment, measured from its left endpoint, by

    q(s)=q₀+s v₀+s² F₀/2,  v(s)=v₀+s[F₀+F(q(s))]/2.

Here q comprises the classical position and the two real left Jacobi positions, with coupled acceleration F=(f(x),f'(x)a,f'(x)b), f(x)=3x exp(−x²/2). Both interpolants hit the exact Verlet nodes. Along a segment
|q'−v|≤D² sup|(F(q))'|/2 and |v'−F(q)|≤D² sup|(F(q))''|/4, componentwise, by Taylor's integral formula. The second identity uses F₀−F(q(s))+s(F(q))'(s). To enclose the derivatives, gather all coarse slabs containing possible segment endpoints; the quadratic position differs from its linear interpolant by at most D²|F₀|/8, its slope is v₀+sF₀ and its second derivative is F₀. The script evaluates the chain rules for f and f'y on these boxes. With the resulting defects ε_q,ε_v, a slab of length t has increments enclosed by

    Δq ∈ t v_start + (t²/2)(F(box) ± ε_v) ± t ε_q,
    Δv ∈ t F(q_start) + (t²/2)(F(q))'(box) ± t ε_v.

These formulas integrate along the interpolant and also hold across fine-node boundaries. The enclosure includes the smooth limit (zero defect), or equivalently take the classical Verlet convergence limit. The right Jacobi solution at a node is (v_b(T)+ib(T))a−(v_a(T)+ia(T))b by unit determinant of the transfer. For the two ordered Green prefix sums, represent trapezoid weights by assigning each physical time to its nearest fine node. Their mass is exact, and the ordered double integral reproduces the diagonal and twice the lower triangle, including time pairs in the same cell. Every possible nearest node is included through the neighboring slab hulls. For a single slab's two-time integral use the rectangle product times its ordered area t²/2. This avoids accumulating a node-count penalty. Positivity of the first-slab b' and subsequent positive b enclosures selects the quadratic branch.

**Action and target phase.** Use carriers exp(ip(x−q)/h) to compare the final coherent states and put the classical action phase on the target ket. Let q_lin denote the piecewise straight classical path through the discrete nodes and S its discrete kinetic-minus-potential action. The continuous reference extremal q_* satisfies the same initial data. First variation with a free terminal position gives the terminal term p_*(2)(q_N−q_*). After subtracting it, the action discrepancy is bounded by the integral of

    (|q_lin'−p_*|² + 3|q_lin−q_*|²)/2 + (D²/4)|q_lin'|².

The last term is the real trapezoidal error for V along a straight segment, using |V''|≤3. For the other terms use the family widths for the enhanced q,v at the *same partial slab time*, plus D²|F₀|/8 for position and D|F₀|/2+D² sup|(F(q))'|/2 for the slope difference. Apply the same estimate separately to the 4096 mesh. Finally replace p_* by (p_N+p_ref)/2 using the enclosed terminal p width times |q_N−q_ref|. This cancellation matters: dividing two raw classical action enclosures by h loses accuracy needlessly.

For two centered normalized packets with their action phases the norm distance is at most √[(Δq²+Δp²)/(2h)] + |ΔS−(p_N+p_ref)Δq/2|/h, by the exact Gaussian overlap. Combine this with the two quantum-to-continuum norm bounds. In the scalar comparison add |ℓ_ref(1+h a_ref)−ℓ_N(1+h a_N)| and |ℓ_ref|·1.6h², then divide by the lower endpoint for |ℓ_N|. This proves the last row of the coverage table.

## Reproduction and verification status

Source checkpoint `535995a912b1479b9222894239151009c7bda79f`. Run with assertions enabled in `experiments/2026-11-24-nova/`, CPython 3.12.14 and NumPy 2.3.5 on IEEE binary64. No interval library beyond the repository's `uniform_049.py` arithmetic is required.

```
python hermite_corrected.py
python direct_corrected.py 1 256
python direct_corrected.py 257 512
python spectral_corrected.py 1024 4096
python bridge_corrected.py
python family_corrected.py 1024 8192
```

Executed both direct tranches, the reference and independent N=1024 checks, the bridge, the infinite family and the Hermite estimate. The direct sweep was also executed in the equivalent pre-refactor layout with the same spectral and scalar operations. The interval parts enclose the classical and polynomial quantities. For the complex amplitude the certificate uses floating-point butterflies with the proved normwise allowance, domain and alias bounds, and the mesh-comparison theorem. Refinement agreement is not used as an error bound. The conventional proofs and arithmetic budgets await independent review; no Lean run or Lean theorem for this result. No displaced-target, certified FFT speedup or continuum-amplitude result is claimed.
