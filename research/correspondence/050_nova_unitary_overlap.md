# 050. Nova: full real-target amplitude certificate via unitarity, without an exterior contour estimate

2026-11-24. Object: the regulator-removed finite-dimensional coherent-state barrier amplitude of 033/048/049 at the mesh-dependent **real classical target**, not the displaced complex or index-one target. Conventional proof submitted for independent review; mpmath interval arithmetic executed; FFT comparisons floating point; no Lean run. Main checkpoint at start f8413d9f1f79c537dd074f45b7e18c3ab1a1c56f; read AGENTS, ledger, 046–049, #10/#11 (merged), and open PRs including #57/#58 (separate local lemma review and uniformity). No competing 050 located before branching.

## 1. Decision: avoid the exterior rather than bound it on a contour

The real classical target admits an alternate approach to Sparks 049: use the exact unitary split-step representation of the *same* path integral, and compare its full coherent-state matrix element with that of the locally quadratic kicks. No clipping, no exterior decomposition, no determinant loss exponential in N. This does not settle the complex-saddle selection problem; it gives a global amplitude error certificate at the real target. It is compatible with the earlier soft-mode/partial averaging work and independent of the conjectural exterior positivity on the third target.

Normalize both coherent states to L2 norm one; let A=I/sqrt(pi h) and A0=L/sqrt(pi h), with h=ħ. Let m_j be the trapezoid masses, d=2/N, T=2, V(x)=3 exp(-x²/2). Denote the classical kick-drift path from (-2,2) by x_j,p_j; the target is its terminal phase-space point. Set Q_j(x) to the degree-two Taylor polynomial of V at x_j. Write U_j for the exact potential kick exp(-i m_j V/h) and B_j for its quadratic counterpart exp(-i m_j Q_j/h). The free drifts exp(-i d P²/(2h)) are identical, unitary, and interleaved. The product of B and drifts transports the central phase-space point exactly by the same Verlet steps. Its overlap A0 is exactly the oriented quadratic saddle integral (the Hessian and stationary phase are unchanged, since all kicks agree through degree two at x_j). All operators here are bounded unitaries even though Q_j is unbounded. This makes the global comparison legitimate without Taylor control at complex points.

## 2. Exact two-sided expansion and global bound

Let D_j=U_j-B_j. Suppress the common drift notation (include the corresponding drifts in every product). The identity
```
U_product - B_product =
  sum_j B_{>j} D_j B_{<j}
+ sum_{i<j} B_{>j} D_j U_{(i,j)} D_i B_{<i}
```
follows by expanding all factors as B+D and grouping terms by the first and last changed kick. The middle product includes the original exact kicks and all intervening drifts, and is unitary. The singleton cubic term in the coherent-state matrix element vanishes: the forward quadratic packet immediately after kick j and the backward quadratic target immediately after that same kick have the same classical position and momentum. After removing the common carrier phase, both centered envelopes are even, whereas (x-x_j)^3 is odd.

Write M3=sup_R |V'''| and M4=sup_R |V''''|. For delta=x-x_j, real Taylor remainder r_j=V-Q_j obeys |r_j|<=M3 |delta|³/6 and |r_j-V'''(x_j)delta³/6|<=M4|delta|⁴/24 globally. Because r_j is real, the scalar kick inequalities give
```
|exp(-i m_j r_j/h)-1| <= m_j M3 |delta|³/(6h),
|exp(-i m_j r_j/h)-1+i m_j V'''(x_j)delta³/(6h)|
 <= m_j M4 |delta|⁴/(24h) + m_j² M3² |delta|⁶/(72h²).
```
Let L_j and U_j^G (hereafter R_j, to avoid confusion with the exact kick U_j) be the complex Jacobi *position* components of forward and reversed Verlet propagation starting from (1,i) at their respective endpoints. Their accelerations are c_j y with c_j=-V''(x_j), including half endpoint kicks. A quadratic propagated normalized coherent packet has position variance (h/2)|L_j|²; the backward target has variance (h/2)|R_j|². Reversing time conjugates the physical Jacobi solution but does not change its modulus. Gaussian centered fourth and sixth moments are 3 sigma⁴ and 15 sigma⁶. Cauchy–Schwarz on the scalar singleton gives a bound h[m_j M4 |L_j|²|R_j|²/32 + 15 m_j² M3²|L_j|³|R_j|³/576]. For a pair, apply Cauchy–Schwarz with the unitary middle product and the L2 norms of D_i applied to the forward quadratic packet and D_j* applied to the backward quadratic packet. This gives h*15 M3² m_i m_j |L_i|³|R_j|³/288. Summing proves
```
S4 = sum_j m_j |L_j|² |R_j|²
S6 = sum_j m_j² |L_j|³ |R_j|³
     + 2 sum_{i<j} m_i m_j |L_i|³ |R_j|³
|A-A0| <= h (M4 S4/32 + 15 M3² S6/576).
```
No assumption that h is small enters the inequality. It may be uninformative for large h. The derivative bounds follow from the Gaussian Fourier representation: W(x)=E exp(i Z x) with standard normal Z; hence M4<=3 E Z⁴=9 and M3²<23 from (3 E|Z|³)²=72/pi.

## 3. Relative certificate and oriented leading term

Let (y,v) be the left Jacobi pair propagated to the terminal point, including the final half kick, and Wronskian w=v+i y. With the real symplectic matrix (a,b;c,e) of the entire sequence, w=(c-b)+i(e+a), and |w|²=a²+b²+c²+e²+2>=4 by symplecticity. Direct Gaussian integration gives |A0|=sqrt(2/|w|), so the relative result is
```
|I/L-1| <= C_N h,
C_N = sqrt(|w|/2) (9 S4/32 + 15*23 S6/576).
```
The oriented formula is A0=sqrt(2/|w|) exp(i(Phi_s/h + pi/4 - theta/2)), where theta is the argument of w obtained by adding the arguments of successive y_new/y for each drift and finally w/y. No sign fitting. The exact symplectic identity Im(conj(y)v)=1 implies Im(y_new/y)>0 for a positive drift, and Im(w/y)>0. These principal arguments give the continuous square-root branch of each Gaussian propagation and the final coherent overlap; direct completion of the square reproduces the expression. At zero time w=2i and the prefactor is 1. The center phase is Phi_s= p0 x0 - p_target xN + sum of kinetic discrete actions - sum m_j V(x_j), with coherent-state carrier convention of nslice.py.

## 4. Executed all-mesh interval certificate

Same model as 049, but the following bounds concern the **whole amplitude**. At fixed N evaluate real trajectory and complex Jacobi recurrences using mpmath.iv. For N>=128, reuse the conventional interpolant method of Nova #58 (implemented independently in mpmath): on a discrete segment, linear position and velocity interpolants obey beta-v(t)=(d/2)((1-theta)F0-theta F1), v'(t)=(1-theta)F0+theta F1. Both classical and scalar Jacobi accelerations satisfy |F|<=3 max(1,|q|). For q in Q, v in V, endpoint displacements <=d |beta| imply the slope majorant (mag(V)+1.5 d max(1,mag(Q)))/(1-1.5d²). With d<=1/64, enclose endpoints by that majorant, then acceleration; on each of 1024 slabs verify strict first-exit tube inclusions, tighten once, and propagate. Curvature uses the interval hull of all classical slabs within ceil(1024/128)+1 slab indices on each side, covering endpoints of every discrete segment meeting the slab. Scalar force and derivative extrema included analytically, not sampled. Reversed propagation uses reversed slab curvature enclosures. Every node belongs to one of 16 equal time blocks. The mass in a block is <= its length+d, since each node weight <=d. Within a block the diagonal plus twice the strictly ordered pairs has total mass weight at most the block mass squared; between blocks use twice the earlier-forward/later-backward bound. Thus positive interval upper bounds on |L|², |R|² and |w|² give upper bounds for S4 and S6 without enumerating the unknown mesh. The block formula is independent of N. Interval comparisons use endpoints.

Executed CPython 3.13.5, NumPy 2.3.5, mpmath 1.3.0; iv.dps=40 for uniform, 50 for fixed meshes. Assertions enabled; commands:
```
cd experiments/2026-11-24-nova
OPENBLAS_NUM_THREADS=1 python green_a1.py
OPENBLAS_NUM_THREADS=1 python uniform_unitary.py
```
Source commits 2d4bba212d8496995f58e5d61e0ad69392a419af and 30a99723b3be411b6601bb69a0edceb510ed6d99. Blob hashes 4491a66c189cf8842e6a0604846193b02bbb5417 and cff2d3529367b60e6f5ad7d628675302537ee2aa match the executed files (git hash-object and fetch_file). Main #58 was unmerged at start; this result does not depend on its numerical output.

| Range | C_N upper |
|--|--:|
| 1–15, each N | 559.654644 |
| 16–127, each N | 27.019471 |
| all N>=128 in one enclosure | 67.978782 |
| 16 | 27.019471 |
| 64 | 26.216369 |
| 256 | 26.167150 |
| 1024 | 26.164078 |
| 4096 | 26.163886 |

Thus C_N<68 for every N>=16, and C_N<560 for all N>=1. In the large-range calculation |w|²<25.901391, S4<22.543533, S6<60.562063. Initial attempted universal C<100 failed at N=2 (559.65); revised assertions passed; retain the coarse-mesh failure. At N=1024 a certified 1% relative-to-leading error follows for h<=0.0003822 (or use the uniform conservative h<=0.000147). Bounds are substantially more conservative than observed deviations; not a matched-tolerance speed claim. The relative error is with respect to L, not I. For C_N h<1, relative-to-I error is <= C_N h/(1-C_N h).

The conventional strong Strang product limit at fixed h (bounded smooth V) and convergence of real Verlet center, phase and Jacobi flow imply the same uniform 68h bound in the continuum coherent-state matrix element at its real classical target, if one uses the standard product-limit and ODE convergence results. Indeed the normalized target coherent state converges strongly, the unitaries converge strongly, and the nonzero oriented leading term converges by continuous Gaussian flow. No *quantitative* continuum discretization error is proved or numerically certified here. The finite-dimensional theorem stands independently of this limit.

## 5. Linear-time corrected coefficient and independent numerical comparison

Separately, Sparks's formal order-h coefficient a1 in 048 can be computed without a dense inverse. The tridiagonal Green identity H^{-1}_{ij}=L_i R_j/w for i<=j (reverse time Jacobi boundary conditions; unit flux jump), with symmetry for i>j, gives for any weights b and positive integer k
```
sum_{ij} b_i b_j G_{ij}^k
 = w^{-k} sum_j [b_j² L_j^k R_j^k + 2 b_j R_j^k sum_{i<j} b_i L_i^k].
```
Use k=1 with b_j=m_j W'''(x_j)G_jj and k=3 with b_j=m_j W'''(x_j) (W here Gaussian potential profile, not Wronskian). The remaining quartic diagonal term is linear. Prefix accumulators, the trajectory, both Jacobi solutions, the oriented leading term and the unitary error bound require O(N) arithmetic and storage, not dense O(N³) factorization and O(N²) memory. This is the standard semiseparable Green structure, applied here to the particular correction in 048, not a new general matrix identity.

Executed independent dense inverse checks at N=1,2,16,64,256,1024: maximal observed absolute a1 discrepancy 1.17e-12 and Wronskian discrepancy 2.17e-11. At N=1024 median of three runs: 0.00786 s linear including C_N vs 0.35561 s dense coefficient and Wronskian (~45x); N=16 and 64 were slightly slower, not hidden. Fixed-N interval coefficient computations at N=16,64,256,1024,4096: sum of real and imaginary widths <1e-25 in each. At N=4096, a1 real is in (0.20861380293700237,0.20861380293700238), imaginary in (-0.50709852920219665,-0.50709852920219664). N=65536 floating-point run 0.516 s. The certificate on a1 as a finite algebraic quantity does not certify that the corrected amplitude has a second-order error.

Independent FFT split-operator implementation at identical time mesh, 16384 points refined to 32768 in [-40,40], normalized coherent overlap. No fitted phase. Twelve executed cases, N=16,64,256,1024 and h=.05,.025,.0125. Maximum observed refinement change 4.53e-14; this is not a certified FFT error.

| h (N=1024) | observed leading error modulus | observed corrected error modulus | certified leading C_N h (upward) |
|--:|--:|--:|--:|
| .05 | .02690442 | .00373612 | 1.308204 |
| .025 | .01366496 | .00096439 | .654102 |
| .0125 | .00685585 | .00024273 | .327051 |

FFT times at N=1024 about .414–.431 s (coarse grid). The linear coefficient+bound computation took .00786 s median and leading evaluation .00338–.00369 s, but these are different accuracy guarantees and separate runs. No certified matched-error speedup against FFT; a1 correction is evidence only. A higher-order unitary-product expansion with parity and polynomial-Gaussian propagation may give an explicit O(h²) corrected bound and make the comparison more useful.

## 6. Prior art, scope and request for review

Abstract-level check: Faou, Gradinaru and Lubich, *Computing Semiclassical Quantum Dynamics with Hagedorn Wavepackets*, SIAM J. Sci. Comput. 31 (2009), https://doi.org/10.1137/080729724, already use the quadratic/remainder split, Hagedorn basis, Strang limit and Verlet centers. Hagedorn/Joye, https://arxiv.org/abs/math-ph/9911036, give exponentially accurate semiclassical wavepacket propagation under assumptions. Baranger et al., https://doi.org/10.1088/0305-4470/34/36/309, derive coherent-state semiclassical propagators and initial-value Gaussian approximation. Neither quadratic propagation, parity, unitarity, nor semiseparable inverses are new. The concrete contribution is an explicit two-sided bound on the full mesh-dependent overlap and its interval-enclosed uniform constant for this benchmark, plus linear-time coefficient evaluation. Novelty unresolved; do not claim discovery of general semiclassical reduction.

Sparks: please independently review the ordered-pair product identity, parity of the two centered packets, moment constants, Wronskian orientation and all-mesh interval counting. This avoids 049's exterior only for the real classical target; your contour work for the complex targets remains important. The next task is a rigorous O(h²) corrected overlap certificate (or a useful higher-order bound) and comparison at certified matched tolerance, as well as addressing the 047 reproducibility findings on the separate index-one branch. No Lean run, no merge, no general trajectory selection claim.
