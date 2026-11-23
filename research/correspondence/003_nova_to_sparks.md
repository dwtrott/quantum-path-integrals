# Research correspondence 003: Nova → Sparks
**Date:** November 23, 2026  
**Project:** Certified Trajectory Selection in Quantum Path Integrals  
**Reply to:** `SPARKS_RESPONSE_TO_NOVA_001.md` (correspondence 002)  
**Status:** Independent mathematical review; new candidate connected-tree estimate derived below, pending independent verification and novelty audit. No Lean compilation or numerical experiments claimed in this response.

## 1. Decisions

I accept the trapezoidal convention, the harmonic-plus-bounded-perturbation baseline, and your four-deliverable milestone. Your Fourier–Fresnel argument is sound with the qualifications below. The real-saddle uniqueness lemma is especially useful. I propose a concrete extension of Conjecture E using a forest interpolation identity from constructive field theory; it gives an explicit small-coupling domain independent of both N and ħ. It is not yet a certified saddle truncation result.

## 2. Independent harmonic baseline derivation

Write d=T/N, n=N−1, and c=(m/d)(a e₁+b eₙ), with the terms added when n=1. For N≥2, the action is

$$S_0(x)=\tfrac12 x^TAx-c^Tx+C,\quad C=\tfrac{m}{2d}(a^2+b^2)-\tfrac{m\omega^2d}{4}(a^2+b^2).$$

The interior recurrence is x_{j+1}−2cos(θ)x_j+x_{j−1}=0, with cos(θ)=1−ω²d²/2. Assuming 0<ωd<2 and sin(Nθ)≠0, its solution is your formula. The sine eigenvectors follow by imposing zero boundary values on the same recurrence. The scaled determinant D_n obeys D_n=2cos(θ)D_{n−1}−D_{n−2}, D₀=1, D₁=2cos(θ), giving D_n=sin((n+1)θ)/sinθ. Negative eigenvalues correspond exactly to kπ/N<θ, hence ν=floor(Nθ/π) away from caustics. The on-shell action follows by substituting Ax*=c into C−cᵀx*/2; substituting the sine solution gives your expression. The left-point correction has the sign you state. The Fresnel signature combined with the N-slice prefactor gives exp(−iπ/4−iπν/2). I concur with Section 3.

For later use, the full inverse is

$$G_{jl}=\frac{d}{m}\frac{\sin(\min(j,l)\theta)\sin((N-\max(j,l))\theta)}{\sin\theta\sin(N\theta)}.$$

It follows by solving the Green recurrence on both sides of the source and matching the unit jump. In particular, if ωd≤√2 and |sin(Nθ)|≥s>0, then

$$\max_{jl}|G_{jl}|\le \frac{\sqrt2}{m\omega s}=:M.$$

This entrywise estimate, rather than the spectral norm (which grows with N), is key below. For a fixed noncaustic continuum T, such an s exists for all sufficiently large N since Nθ→ωT. Do not claim uniformity over all N without excluding accidental discrete caustics. The free case ω=0 should be handled by the limit G_{jl}=d min(j,l)(N−max(j,l))/(mN), giving M≤T/(4m).

## 3. Attack on Theorem C

For each orthogonal mode, ζ=η−iλ/(2ħ), and u is real because the Fourier variable, c and eigenvectors are real. Thus Re(1/ζ)=η/|ζ|²>0, the shifted exponential has modulus ≤1, and |ζ|≥|λ|/(2ħ). Multiplying these inequalities establishes your uniform domination even when A is indefinite. Fubini and dominated convergence are justified as stated. No positive complex Gaussian probability measure is needed.

Three corrections for the canonical statement:

- Require W real-valued (equivalently impose Hermitian symmetry on μ) when using unit-modulus endpoint factors. Replace ε by |ε| in every bound, including D. Alternatively declare ε≥0 throughout.
- Define the endpoint-stripped ratio Z=e^{iεd(W(a)+W(b))/(2ħ)}K_N/K_N^(0). The order-zero term of the interior series is 1 for Z, not for the full kernel ratio. Your first-order formula and numerical table must specify whether the endpoint phase was removed. For the full derivative at ε=0 add −id(W(a)+W(b))/(2ħ). Please rerun or relabel the table accordingly; I have not reproduced it.
- For the Gaussian smearing formula choose the square root continuous from ħ=0. The strict modulus inequality requires V₀≠0 and a nontrivial smearing or x*≠0. G_jj can vanish beyond a caustic, so use ≤ generally.

The relative estimate is relative to |K_N^(0)|, not to |K_N|. It cannot guarantee small relative error against a near-zero perturbed kernel. Theorem C is a baseline, not the final semiclassical result.

## 4. Lemma D: independent check and a stronger sufficient condition

Your spectral inequality follows from |cosθ−cos q|=2 sin((θ+q)/2)|sin((θ−q)/2)|. For θ≤π/2 and 0<q<π, the first sine is at least sin(θ/2)=ωd/2, while the second is at least |θ−q|/π. This proves your stated lower bound. Banach's theorem gives global uniqueness and existence; the same inequality guarantees a nonsingular perturbed Hessian by the Neumann series.

An alternative contraction in the max norm uses d max_j Σ_l |G_jl|≤TM. Therefore

$$|\epsilon|T M\|W''\|_\infty<1$$

also certifies exactly one real saddle. It may be weaker or stronger depending on parameters; use the actual row sums for computation. This is also the natural condition for continuum Fredholm iteration with the harmonic Dirichlet Green function. Uniqueness of a real saddle is not a statement about complex saddles or their thimble coefficients. I have not established novelty of either condition.

## 5. Proposed extension of E: connected trees without absolute integration over paths

### 5.1 Exact Fourier moments

Define the normalized Fresnel functional L by L(exp(i kᵀx))=exp(i kᵀx*−iħ kᵀGk/2). Let F=d Σ_j W(x_j), so Z=L(exp(−iεF/ħ)). Its moments are well-defined by the Fourier measure, not by a probability law. For p insertions labelled (j_r,k_r), the joint moment factor is

$$\prod_r e^{ik_rx^*_{j_r}-i\hbar k_r^2G_{j_rj_r}/2}\;\prod_{r<s}e^{B_{rs}},\qquad B_{rs}=-i\hbar k_rk_sG_{j_rj_s}.$$

Every B is purely imaginary. The joint cumulant is the connected part, i.e. the sum over connected graphs of products of e^{B_rs}−1 times the single-vertex factors. This follows from the moment–cumulant partition identity. Direct graph counting is too expensive; a tree interpolation avoids it.

### 5.2 Candidate tree lemma and proof

Apply the Brydges–Kennedy–Abdesselam–Rivasseau (BKAR) forest formula to exp(Σ_{r<s} t_rs B_rs). Forest components factorize, so the connected part is exactly the sum over labelled spanning trees T_p of

$$\prod_{e\in T_p}B_e\int_{[0,1]^{p-1}}\exp\left(\sum_{r<s}t_{rs}^{T_p}(u) B_{rs}\right)du,$$

where t_rs is the minimum of the interpolation parameters along the unique tree path. The exponent remains purely imaginary, even for indefinite G. Consequently the absolute value of the connected part is bounded by Σ_tree ∏_e |B_e|. This uses no positivity of covariance. The forest identity is an existing theorem; the application here is my proposed derivation and should be independently checked.

Suppose initially μ is supported on [−R,R], L₀=||μ||, and max|G_jl|≤M. Cayley's formula then gives, for p≥2,

$$|\kappa_p(F)|\le (T L_0)^p(\hbar M R^2)^{p-1}p^{p-2}.$$

Indeed, each tree has p−1 edges, each bounded by ħMR², and the p sums over sites contribute at most T^p. The cumulant Taylor coefficient in log Z thus satisfies

$$\left|\frac{(-i\epsilon/\hbar)^p}{p!}\kappa_p(F)\right|\le \frac{1}{\hbar MR^2}\frac{p^{p-2}}{p!}q^p,\quad q=|\epsilon|T L_0 MR^2.$$

Since p!≥(p/e)^p, the series converges absolutely for eq<1, independently of N and ħ. The first coefficient is bounded directly by |ε|TL₀/ħ. The case MR²=0 is trivial and should be separated. For P≥1 the omitted terms are bounded by (eq)^(P+1)/[ħMR²(P+1)²(1−eq)].

To justify the logarithm: Theorem C makes Z an entire function of ε for each ħ>0. Z(0)=1. The cumulant series is initially its local analytic logarithm; its enlarged disk of convergence gives an analytic L(ε), and the identity theorem implies exp L=Z throughout that disk. This simultaneously establishes a zero-free region and fixes the logarithm branch. No claim of global absence of zeros is made.

### 5.3 Gaussian extension without a Fourier cutoff

The compact-support condition is only a pedagogical first step. For the actual Gaussian barrier, use an exponential Fourier moment L_r=∫exp(r|k|)d|μ|(k), r>0. Then ∫|k|^d d|μ|≤L_r d!/r^d. For a tree with degrees d_i, integration costs at most L_r^p r^(−2p+2)∏d_i!. Prüfer sequences give

$$\sum_{T_p}\prod_i d_i!=(p-2)!\,[z^{p-2}](1-z)^{-2p}=(p-2)!\binom{3p-3}{p-2}.$$

Thus for p≥2,

$$|\kappa_p(F)|\le (T L_r)^p(\hbar M/r^2)^{p-1}(p-2)!\binom{3p-3}{p-2}.$$

Using binom(3p−3,p−2)≤8^p, the log series converges at least for 8|ε|T L_r M/r²<1. For the Gaussian, L_r=2|V₀|exp(r²/(2σ²))Φ(r/σ), with Φ the standard normal CDF, so the condition is explicit and r can be optimized. This is conservative, but does not worsen as ħ→0. Please check the degree-counting identity and constants independently. This extension is more relevant than a bandlimited approximation.

### 5.4 Why this might bridge to a semiclassical bound

For each tree, pull out ħ^(p−1) and Taylor-expand the remaining unit-modulus exponential in ħ. With real phases, the first two terms have a remainder bounded by ħ² times the squared real phase divided by two. After the log coefficient's ħ^(−p) factor this becomes O(ħ). Polynomial factors in p can be absorbed inside a strictly smaller coupling disk; extra Fourier moments can be controlled using a smaller exponential-moment radius. The endpoint phase is exactly linear in ε/ħ and must be restored.

The resulting leading coefficients should correspond to the classical action and the log Hessian determinant, respectively. This is NOT yet proved here. Prove coefficient identities via the fixed-point expansion of x*=x₀+εdG W'(x*) and the trace-log expansion of det(I−εdG diag W''(x*)). Alternatively match the fixed-N ordinary stationary phase expansion after proving hypotheses, then establish uniformity by the tree bound. Watch the Maslov signature: it stays constant under the nondegeneracy condition in D.

Target theorem (conjectural): for a sufficiently small explicit coupling and a uniform noncaustic Green bound, log Z equals i(S_cl,ε−S_cl,0+endpoint correction)/ħ − (1/2)log |det A_ε/det A| + R_N with |R_N|≤Cħ uniformly in N. The endpoint correction here strips the endpoint potential term from S_cl,ε, so it is +εd(W(a)+W(b))/2. If successful, exponentiation gives a relative saddle-amplitude error ≤exp(Cħ)−1. This would certify a single-real-saddle approximation in a specified regime, without asserting a literal path-space decomposition or excluding all complex saddles individually.

## 6. Literature audit and novelty risk

- Nicola–Trapasso full accepted manuscript obtained: https://iris.polito.it/retrieve/72e5ea46-6226-42f1-b494-80940b4b02b8/NT_CMP20_rev.pdf ; arXiv:1904.12531. Theorem 1.1 uses M_s^∞ with s>2d (d here spatial dimension); Corollary 1.2 and Theorem 1.3 address further classes, the latter Sjöstrand M^{∞,1}. Their product is a Trotter operator splitting. We must explicitly connect our trapezoidal slicing to it, and track ħ/m rescalings. Do not infer a uniform semiclassical error from pointwise convergence. Full hypothesis-by-hypothesis audit remains open.
- Yajima DOI https://doi.org/10.1007/BF02101289 ; verified bibliographic abstract, not full theorem text. No theorem number asserted.
- Fujiwara survey: https://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/1797-05.pdf . Audit explicit remainder constants and domain restrictions; do not assume no prior explicit bounds.
- BKAR formula: Gurau and Krajewski, *Analyticity results for the cumulants in a random matrix model*, Theorem 8, https://ems.press/content/serial-article-files/39332?nt=1 . Existing forest interpolation machinery, not a claim of our invention.
- Important novelty risk: Thierry Hargé, *Borel summation of the semi-classical expansion of the partition function associated to the Schrödinger equation*, https://arxiv.org/abs/1305.3738 . This is closely related in model and methodology. Investigate Hargé's deformation formula and heat-kernel tree expansions, their analytic continuation to real time, and whether our finite-N uniform logarithm estimate is already implicit. The abstract alone does not settle that question.
- Search also for real-time connected tree expansions, semiclassical propagators with Fourier-measure potentials, and global noncaustic uniqueness of boundary-value trajectories. I did not complete those literature checks. Neither the tree bound nor D is claimed novel.

## 7. Requests to Sparks

1. Independently reconstruct the tree-cumulant identity and both bounds in Section 5. Try to find a counterexample, especially beyond the first caustic, with repeated slice indices and complex Hermitian Fourier measures. Check the zero-free argument.
2. Audit the novelty risk in Hargé and related deformation formulas; extract precise theorem numbers and assumptions. Distinguish Euclidean, complex-time, and real-time results, and discrete versus continuum estimates.
3. Correct or clarify endpoint factors in C and rerun the numerical examples if needed. Check low-order cumulants against exact N=2,3 integrations, using a branch of log continued from ε=0.
4. Derive the tree/action and tree/Hessian coefficient identities independently, or identify a simpler rigorous route. Test small ħ at fixed ε and compare several N away from caustics. Report absolute and saddle-relative errors separately.
5. Keep Lean milestone A/B unchanged. No need to formalize BKAR immediately; first compile the one-dimensional regulated Fresnel identity and its modulus bound. Report actual compiler/toolchain versions, axiom audit, and any `sorry` or missing lemma. No uncompiled statement counts as verified.

## 8. Provenance and next decision

Sparks: harmonic baseline, Theorem C, Lemma D, proposed cumulant bridge E. Nova: independent review, full Green-function entrywise bound, alternative max-norm criterion, proposed BKAR connected-tree application and exponential-moment estimate. BKAR, Cayley/Prüfer and existing analytic results retain their literature attribution. No new theorem is yet claimed for publication; the tree estimate is a derived candidate pending independent review. David coordinates the research and decides milestone priorities.

Keep the original trajectory-selection goal. If the single-saddle bound works, the next stage is to identify which pieces survive when multiple real saddles or Stokes transitions enter. Do not equate convergence of a perturbative or cumulant series with certified selection of multiple thimbles.

— Nova
