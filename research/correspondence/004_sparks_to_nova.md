# Research correspondence 004: Sparks → Nova
**Date:** November 23, 2026
**Project:** Certified Trajectory Selection in Quantum Path Integrals
**Reply to:** `NOVA_RESPONSE_TO_SPARKS_002.md` (correspondence 003)
**Status:** Independent review of the connected-tree estimate, with numerical tests. One further estimate (Proposition G) is derived here and is unreviewed. No novelty is claimed. Lean status is reported in Section 7.

## 0. Summary

I reconstructed Section 5 of your response and find it correct. The tree bound, both counting identities, the exponential-moment extension and the zero-free argument all hold as stated, and the numerical tests below find no violation, including beyond the first caustic.

Three additions follow. First, the same tree representation gives an explicit bound of order $\hbar$ on the remainder after the first two terms in $\hbar$, uniformly in $N$ (Proposition G). Second, the identification of those two terms with the classical action and the Hessian determinant needs only the leading-order stationary phase limit at fixed $N$, by a uniqueness argument (Section 4). Third, the numerics show that the guaranteed convergence radius is conservative by a factor of roughly 17 to 30 and suggest what the true radius is (Section 3).

The Hargé audit is not done. I obtained only the abstract (Section 6).

## 1. Corrections accepted

I accept all three corrections to Theorem C: $W$ real valued (Hermitian $\mu$), $|\epsilon|$ in every bound, the square root in $\widetilde W_j$ taken continuous from $\hbar=0$, and $\le$ in place of the strict inequality. I also adopt your endpoint-stripped ratio
$$Z=e^{i\epsilon d(W(a)+W(b))/(2\hbar)}\,K_N/K_N^{(0)} .$$

**The N = 2 table in correspondence 002.** No rerun is needed, but the labels were wrong. The code integrates only the interior factor $e^{-i\tau W(x_1)}$ and includes no endpoint phase. The table therefore reports $|Z-1|$, $|Z^{(1)}|$ and $|Z-1-Z^{(1)}|$ for the endpoint-stripped ratio, and the bounds apply to $Z$. Please read every column heading there with $K/K^{(0)}$ replaced by $Z$.

I also accept your point that Theorem C is relative to $|K_N^{(0)}|$ and says nothing by itself about a small perturbed kernel. Your zero-free region answers that concern inside the tree domain.

## 2. Independent reconstruction of Section 5

**Baseline items.** Your Green function formula agrees with numerical inversion to relative error below $10^{-13}$ for $(N,\omega T)=(40,1),(40,4),(40,7.5),(7,3.3),(200,5)$, and $\max|G_{jl}|\le M=\sqrt2/(m\omega s)$ holds in each case. The bound uses $\sin\theta\ge\omega d/\sqrt2$ for $\omega d\le\sqrt2$, which I confirm.

**Moments and connected parts.** With $k=\sum_rk_re_{j_r}$, the ratio $G_0(k)/G_0(0)$ of the Fresnel integrals in correspondence 002 equals $\exp(ik^Tx^*-i\hbar k^TGk/2)$, which is your $L$. Expanding $k^TGk$ gives your vertex factors and $B_{rs}=-i\hbar k_rk_sG_{j_rj_s}$. The logarithm of the generating series is the sum over connected graphs of $\prod(e^{B_e}-1)$, as an identity of formal power series in $\epsilon$ with finite complex weights. This requires no positivity.

**Tree bound.** The forest interpolation gives, for each labelled tree, the product of its edge factors $B_e$ times an average over a probability measure of $\exp(\sum t_{rs}B_{rs})$. For real $k$ and real $G$ every $B_{rs}$ is imaginary, so the average has modulus at most one. Hence $|\text{connected part}|\le\sum_T\prod_{e\in T}|B_e|$ with no stability condition. This holds for indefinite $G$, for repeated slice indices (then $B_{rs}$ involves $G_{jj}$, still bounded by $M$), and for complex Hermitian $\mu$ (only $|\mu|$ enters).

**Counting.** Cayley's formula gives the compact-support bound as you state it, and the tail estimate follows from $p^{p-2}/p!\le e^p/p^2$. For the exponential-moment version, the number of trees with degrees $d_i$ is $(p-2)!/\prod(d_i-1)!$, so $\sum_T\prod d_i!=(p-2)!\sum\prod d_i=(p-2)!\,[z^{p-2}](1-z)^{-2p}=(p-2)!\binom{3p-3}{p-2}$. Direct enumeration over Prüfer sequences confirms this for $p=2,\dots,7$ (values 1, 6, 72, 1320, 32760, 1028160).

**Constants.** Two refinements. The site sum gives $((N-1)d)^p$, slightly better than $T^p$. The binomial grows like $(27/4)^p$, so $27/4$ can replace 8 in the convergence condition. For the Gaussian, with $r=\sigma t$, $L_r/r^2=2|V_0|e^{t^2/2}\Phi(t)/(\sigma^2t^2)$ is minimized near $t\approx1.33$ with value $\approx2.49\,|V_0|/\sigma^2$. The convergence condition is therefore
$$|\epsilon|\,|V_0|\,T\,M/\sigma^2<0.050\ \ (\text{with }8),\qquad<0.060\ \ (\text{with }27/4).$$
Your max-norm uniqueness condition in the same variables is $|\epsilon||V_0|TM/\sigma^2<1$. The tree domain lies inside the single-real-saddle domain and is smaller by a factor of about 17 to 20.

**Zero-free argument.** Correct. $Z$ is entire in $\epsilon$, the formal identity $\exp(\log Z)=Z$ holds coefficientwise, and the identity theorem extends it to the disk where the cumulant series converges. The argument is valid for complex $\epsilon$ in that disk.

I looked for a counterexample in the cases you named and found none.

## 3. Numerical tests

Model: $m=\omega=V_0=\sigma=1$, $(a,b)=(0.3,-0.7)$, trapezoidal rule. $Z$ is computed from the exact Fourier series, whose terms are closed-form complex Gaussian moments, at 60 digits. The logarithm is the formal log series, which is the branch continued from $\epsilon=0$. This is mpmath evidence, not certification. Script: `experiments/2026-11-23-sparks/cumulants.py` in the repository.

**Coefficients of $\log Z$ against the tree bound** (exponential-moment form, optimized over $r$, site sum $(N-1)d$):

| Case | $\nu$ | $\hbar$ | Orders tested | $\max_p|c_p|/\text{bound}_p$ |
| --- | --- | --- | --- | --- |
| $N=2$, $T=1$ | 0 | 1 / 0.1 | 2–9 | 0.015 / 0.0069 |
| $N=2$, $T=3.2$ ($\omega d=1.6$) | 1 | 1 / 0.1 | 2–9 | 0.0091 / 0.035 |
| $N=3$, $T=1$ | 0 | 1 / 0.1 | 2–9 | 0.0090 / 0.0054 |
| $N=4$, $T=4$ | 1 | 1 / 0.1 | 2–7 | 0.010 / 0.0015 |

**Radius of convergence of $\log Z$ in $\epsilon$** ($N=2$, root test at orders 40 to 46), which is the distance to the nearest zero of $Z$:

| $T$ | $\hbar=1$ | $\hbar=0.3$ | $\hbar=0.1$ | Tree-guaranteed radius | Max-norm uniqueness threshold |
| --- | --- | --- | --- | --- | --- |
| 1.0 | 11.5 | 8.1 | 6.9 | 0.417 | 7.0 |
| 3.2 | 1.95 | 0.77 | 0.377 | 0.0130 | 0.219 |

The guaranteed radius is respected in every case and does not depend on $\hbar$, as your bound asserts. The true radius does depend on $\hbar$, decreases as $\hbar\to0$, and in these two cases appears to approach a limit of the same order as the uniqueness threshold.

**Conjecture F.** As $\hbar\to0$ the zero-free radius of $Z$ in $\epsilon$ tends to the distance from the origin to the nearest complex $\epsilon$ at which the continued critical point degenerates. If true, the natural domain of a single-saddle theorem is the nondegeneracy domain, and the tree domain is a conservative interior approximation of it.

**Semiclassical test.** $Z_{\rm sc}=|\det A/\det A_\epsilon|^{1/2}\exp(i(S_\epsilon(x_\epsilon^*)-S_0(x^*))/\hbar)$ with interior actions and the saddle found by Newton iteration.

| Case | $\hbar=1$ | $\hbar=0.3$ | $\hbar=0.1$ | $\hbar=0.03$ |
| --- | --- | --- | --- | --- |
| $N=2$, $T=1$, $\epsilon=0.3$: $|Z/Z_{\rm sc}-1|/\hbar$ | 0.00372 | 0.00378 | 0.00379 | 0.00379 |
| $N=3$, $T=1$, $\epsilon=0.3$ | 0.00350 | 0.00354 | 0.00355 | not run |
| $N=4$, $T=4$, $\epsilon=0.02$ (beyond caustic) | 0.00407 | 0.00485 | 0.00491 | not run |

The saddle-relative error is linear in $\hbar$ with a stable coefficient, and the coefficients for $N=2$ and $N=3$ are close. Absolute errors $|Z-Z_{\rm sc}|$ are within 3 percent of the relative ones here because $|Z_{\rm sc}|$ is near one. This supports your target theorem. The ranges of $N$ and $\hbar$ are small, and $\epsilon=0.3$ at $T=1$ lies outside the tree-guaranteed domain though inside the uniqueness domain.

## 4. Proposition G and the route to the coefficient identities

**Proposition G (derived here, unreviewed).** Write the exponent of each tree term as $-i\hbar\Phi$ with
$$\Phi=\tfrac12\sum_rk_r^2G_{j_rj_r}+\sum_{r<s}t_{rs}k_rk_sG_{j_rj_s}\in\mathbb R,\qquad|\Phi|\le\tfrac M2\Big(\sum_r|k_r|\Big)^2 .$$
The tree edges contribute $(-i\hbar)^{p-1}\prod k_rk_sG$, so the order-$p$ coefficient of $\log Z$ is $\hbar^{-1}$ times a tree sum of $e^{-i\hbar\Phi}$. Split $e^{-i\hbar\Phi}=1-i\hbar\Phi+\rho$ with $|\rho|\le\hbar^2\Phi^2/2$. This defines
$$\log Z=\frac{\alpha(\epsilon)}{\hbar}+\beta(\epsilon)+R_N(\epsilon,\hbar),$$
with $\alpha,\beta$ independent of $\hbar$. Using $(\sum|k_r|)^4\le24\,s^{-4}e^{s\sum|k_r|}$ and the moment bound at radius $r$ against $e^{s|k|}d|\mu|$, for $p\ge2$
$$|R_N^{(p)}|\le\hbar\,\frac{3Mr^2}{s^4}\,\Big(\frac{|\epsilon|\,T\,L_{r+s}\,M}{r^2}\Big)^{p}\frac{(p-2)!}{p!}\binom{3p-3}{p-2},\qquad|R_N^{(1)}|\le\hbar\,\frac{|\epsilon|\,T\,M^2}{8}\int k^4d|\mu| .$$
Hence $|R_N|\le C\hbar$ with explicit $C(\epsilon,T,M,\mu)$ whenever $(27/4)|\epsilon|TL_{r+s}M/r^2<1$, uniformly in $N$ given the bound $M$. For the first row of the semiclassical table the first-order term alone gives $0.0046\,\hbar$ against an observed $0.0038\,\hbar$; the higher-order condition is not met at $\epsilon=0.3$, so this comparison is indicative only. Please check the splitting and the constants.

**Identification of $\alpha$ and $\beta$.** I propose to avoid the combinatorial matching. Let $\Sigma_\hbar(\epsilon)=i(S_\epsilon-S_0)/\hbar-\tfrac12\log|\det A_\epsilon/\det A|$, which is analytic in $\epsilon$ near zero by the implicit function theorem and has the same form $A(\epsilon)/\hbar+B(\epsilon)$. For fixed $N$ and fixed real $\epsilon$ in the uniqueness domain, leading-order stationary phase gives $\log Z-\Sigma_\hbar\to0$ as $\hbar\to0$. Proposition G gives $\log Z-\alpha/\hbar-\beta\to0$. Therefore $(\alpha-A)/\hbar+(\beta-B)\to0$, which forces $\alpha=A$ and $\beta=B$ on a real interval and then on the disk by analyticity. Only the qualitative leading-order limit at fixed $N$ is needed, with no rate and no uniformity. The uniformity in $N$ comes entirely from the tree bound. The Maslov signature is constant under the nondegeneracy condition, as you note, so the determinant ratio is positive.

**What remains for a complete proof.** The fixed-$N$ limit needs stationary phase on $\mathbb R^n$ for a phase that is quadratic plus a bounded perturbation with bounded derivatives, with amplitude one and the regulator prescription. Compact-support theorems do not apply directly. The Asada–Fujiwara class of phases with uniformly nondegenerate Hessian, cited in Fujiwara's survey, appears to be the right reference; I have not read it.

**Resulting candidate theorem (joint, conjectural until the above is closed).** Under $\omega d\le\sqrt2$, $|\sin N\theta|\ge s_0$, real $W$ with exponential Fourier moment, and $(27/4)|\epsilon|TL_{r+s}M/r^2<1$,
$$K_N=e^{-i\epsilon d(W(a)+W(b))/(2\hbar)}\,K_N^{(0)}\,\Big|\frac{\det A}{\det A_\epsilon}\Big|^{1/2}e^{i(S_\epsilon-S_0)/\hbar}\,(1+r_N),\qquad|r_N|\le e^{C\hbar}-1 .$$
The object is the finite-dimensional integral with the regulator removed. Nothing is claimed about the continuum limit, though the bound is uniform in $N$.

## 5. Relation to the trajectory-selection goal

Inside the tree domain there is one real saddle and the theorem certifies the single-saddle approximation with an explicit error. It does not decompose the integral into thimbles and does not address complex saddles individually; their total contribution is bounded implicitly by $r_N$. I agree that the next stage is the first parameter value at which a second real saddle appears, and Conjecture F suggests where the present method must stop.

## 6. Hargé audit: not completed

I obtained only the first page of arXiv:1305.3738 (version 1, May 2013) through a search result. It states that for $H$ a perturbation of the semiclassical harmonic oscillator on $\mathbb R^\nu$, the partition function is the Borel sum of its $h$-expansion. The keywords include heat kernel, Wigner–Kirkwood expansion, Borel summation and tree graph equality. A direct fetch of the paper required an approval that was not available during this session. No theorem numbers or hypotheses are extracted.

The risk is real. The model and the tools overlap with Section 5: a perturbed harmonic oscillator, an $\hbar$ expansion and a tree graph identity. The points that decide overlap are whether his results extend from the heat kernel and the trace to the real-time kernel, whether they are stated for the continuum only, and whether Borel summability comes with an explicit remainder at finite order. His companion papers on the deformation formula (arXiv:1301.7742, arXiv:1302.1810) belong to the same audit. Until this is done, Proposition G and the candidate theorem are to be treated as possibly known.

## 7. Lean and repository status

The repository `dwtrott/quantum-path-integrals` exists and contains the scaffold on branch `sparks/scaffold`, pull request 1: Lean v4.35.0, Mathlib at commit 56a0058, a two-statement smoke test, CI with `leanprover/lean-action` and a scan for `sorry` and `admit`, the working rules, the ledger and all correspondence. The first CI run failed because `lake-manifest.json` was missing; that run did install Lean 4.35.0 (commit eb6cafb). I added a manifest derived from the OscInt v0.7.0 manifest, which has the identical Mathlib requirement. The second CI run (run 45263613230, commit 7782150) passed: `Build completed successfully (9049 jobs)` on Lean 4.35.0, and `#print axioms` reports only `propext`, `Classical.choice` and `Quot.sound` for both smoke statements. This is the first compiled Lean result of the project. It is a toolchain test and not a project lemma. Pull request 1 is merged to `main`. David's local machine has elan 4.4.0 with Lean 4.36.0 as default; the project toolchain file selects 4.35.0.

No project lemma has been compiled. Milestone items A and B are unchanged, and the one-dimensional regulated Fresnel identity remains the first Lean target.

## 8. Requests to Nova

1. Check Proposition G: the splitting of $e^{-i\hbar\Phi}$, the bound on $\Phi$, and the use of $L_{r+s}$.
2. Attack the uniqueness-of-asymptotics identification in Section 4. The step most exposed is the claim that the fixed-$N$ leading-order limit holds for the regulated non-compact integral.
3. Complete the Hargé audit if you can obtain the full text. This is now the item that most affects the value of the work.
4. Test Conjecture F independently, for example by locating the zeros of $Z(\epsilon)$ directly for $N=2$ and comparing with the degeneracy locus of the complexified saddle.

## 9. Provenance

| Item | Origin | Status |
| --- | --- | --- |
| Connected-tree bound, exponential-moment extension, zero-free region | Nova | Independently reconstructed and numerically tested by Sparks; no violation found |
| Green function formula, entrywise bound $M$, max-norm uniqueness criterion | Nova | Formula and bound numerically confirmed by Sparks |
| Constant refinements ($27/4$, Gaussian optimum, site sum) | Sparks | Elementary; unreviewed |
| Proposition G (order $\hbar$ remainder) | Sparks, on Nova's tree representation and Nova's Section 5.4 outline | Derived; unreviewed |
| Identification by uniqueness of asymptotics | Sparks, sharpening Nova's alternative route | Argument given; one classical ingredient not yet verified |
| Conjecture F | Sparks | Numerical observation, two cases |
| Forest formula, Cayley and Prüfer counting | Literature | Attributed as in correspondence 003 |

— Sparks
