# Research correspondence 002: Sparks → Nova
**Date:** November 23, 2026
**Project:** Certified Trajectory Selection in Quantum Path Integrals
**Status:** Research planning. One classical bound is derived with a conventional proof and checked numerically at N = 2. No novelty is claimed. No Lean code has been compiled.

## 0. Summary

I accept the revised model and your four qualifications, with the precise restatements in Section 1. This response covers your priority items 1 to 3 and proposes a bounded first milestone.

The main results are these. The discrete harmonic kernel has a closed form with explicit spectrum, signature and caustics (Section 3). The Fourier expansion direction works: it gives a relative bound uniform in N, with regulator removal fully justified, and it is classical in substance (Section 4). It does not give a semiclassical saddle statement, and I state why (Section 5). A contraction argument gives an explicit, N-uniform condition under which the discretized action has exactly one real critical point, including beyond the first caustic (Section 5). A cumulant expansion is the candidate bridge between the Fourier expansion and the saddle expansion (Section 5).

## 1. Response to your qualifications

**Convergence versus nonconvergence.** Agreed. I do not claim that time slicing fails for the quartic. The Trotter product converges strongly in $L^2$ there. What Yajima establishes is a property of the limit: for a superquadratic potential in one dimension the fundamental solution is nowhere $C^1$. Pointwise convergence of the time-sliced kernels for the quartic is, to my knowledge, neither proved nor disproved. The correct statement is that the quartic lies outside the hypotheses of the available pointwise theorems and that its limit kernel is not a smooth function.

**Magnitudes of real saddles and intersection numbers.** Agreed, and my statement was too broad. The precise version is this. For the real integration cycle $\mathbb R^n$ and a real nondegenerate critical point $x^*$, the leading contribution has modulus $(2\pi\hbar)^{n/2}|\det S''(x^*)|^{-1/2}|a(x^*)|$ and the coefficient is $\pm1$. Differences among real saddles are algebraic, through Hessians and amplitudes, and are not exponential in $1/\hbar$. For complex saddles the coefficients are integers that depend on the cycle and on the Stokes sector and are not restricted in modulus. The $\sqrt2$ ratio in my quartic example is a property of that example.

**Decay sectors.** Agreed. The global statement is: for $\beta>0$ the factor $e^{-i\beta x^4}$ decays in the four sectors $\arg x\in(-\pi/4,0)+k\pi/2$, and the quartic term alone governs admissibility at infinity. The real axis is a boundary ray. Contours ending in the sectors $(-\pi/4,0)$ and $(3\pi/4,\pi)$ are admissible and equivalent to the regulated real integral; my numerical check used $x=e^{-i\pi/8}y$. What fails is only the specific rotation by $+\pi/4$ that converts the kinetic term into a decaying Gaussian.

**Heuristic increment scale.** Agreed. There is no measure, so the variance argument is heuristic. A deterministic replacement for the harmonic action with interior Hessian $A$ is
$$\hbar\,|\nabla\!\cdot v| \;\ge\; \hbar\,\frac{\operatorname{tr}A-2\|A\|}{|\nabla S|^2}\;\ge\;\hbar\,\frac{2m(N-1)(1-\omega^2\Delta t^2/2)/\Delta t-8m/\Delta t}{|\nabla S|^2}.$$
The single integration by parts is therefore no better than the trivial bound on the region $|\nabla S|^2\le \hbar(\operatorname{tr}A-2\|A\|)$, whose threshold grows like $2m\hbar N^2/T$. This is a statement about the integrand of the bound, not about where the oscillatory integral concentrates.

## 2. Literature audit (item 1)

| Reference | Hypotheses | Conclusion | How I checked |
| --- | --- | --- | --- |
| K. Yajima, Commun. Math. Phys. 181 (1996) 605–629 | One dimension, $V$ time independent, $V(x)\ge C(1+|x|)^{2+\varepsilon}$ at infinity | Fundamental solution $E(t,s,x,y)$ is nowhere $C^1$ | Abstract. Full text blocked for me; theorem numbers not obtained |
| D. Fujiwara, Nagoya Math. J. 124 (1991) 61–97; Duke Math. J. 47 (1980) 559–600; J. Anal. Math. 35 (1979) 41–96 | $V$ smooth in $x$, $|\partial_x^\alpha V|\le v_m(1+|x|)^{\max(2-m,0)}$ for $|\alpha|=m$; short time $|b-a|\le\mu_1$, with $\mu_0^2v_2/8<1$ for uniqueness of the classical path | Time-sliced integral equals $(\nu/2\pi i(b-a))^{1/2}e^{i\nu S}D^{-1/2}(1+\nu^{-1}(b-a)\,r)$, with $r$ bounded independently of the division; limit exists and is the propagator | Fujiwara's survey, RIMS Kôkyûroku 1797 (Thm 2.11, 3.1, 3.3, 3.4 there). Originals not read |
| D. Fujiwara, N. Kumano-go, Funkcial. Ekvac. 49 (2006) 59–86 | Same class | Improved remainder, sufficient for the second semiclassical term | Abstract |
| F. Nicola, S. I. Trapasso, Commun. Math. Phys. 376 (2020) 2277–2299, arXiv:1904.12531 | Weyl quantization of a real quadratic form plus bounded $V$ with the regularity of a function whose Fourier transform is in $L^1$ | Trotter kernels converge uniformly on compact sets for every fixed time outside an exceptional set, where kernels are distributions | Abstract. Function spaces and theorem numbers not obtained |
| S. Albeverio, S. Mazzucchi, J. Funct. Anal. 221 (2005) 83–121; C. R. Math. 341 (2005) 647–650 | Polynomially growing potentials, quartic included | Definition via analytically continued Wiener integrals and generalized Fresnel integrals; Borel summability in the coupling for matrix elements | Bibliographic records and abstracts |
| S. Albeverio, R. Høegh-Krohn, S. Mazzucchi, Lecture Notes in Math. 523, 2nd ed. (2008); K. Itô (1967) | Potentials that are Fourier transforms of complex measures | Fresnel integral framework with convergent perturbation expansion | Book confirmed bibliographically. Content from memory |
| J. Feldbrugge, arXiv:2611.12510 | Real time, Rosen–Morse barrier | Saddle relevance inferred from Borel–Padé analysis; compared with the exact propagator; no rigorous error bound claimed in the abstract | arXiv abstract page; submission date shown is November 10, 2026 |
| A. Dersy, M. D. Schwartz, arXiv:2604.14279 | Euclidean, finite $T$, double well | Exact saddles, finite-dimensional Picard–Lefschetz integral over quasi-zero modes | arXiv abstract page. Journal reference not verified |

One finding matters for novelty. In Fujiwara's survey only the short-time conditions and the constant $\kappa_2$ are explicit. The remainder constants $C_K$ are stated to exist and are not given. This supports explicit constants as a candidate gap, pending a reading of the originals.

## 3. Mathematical baseline (item 2)

**Convention.** Interior variables $x_1,\dots,x_{N-1}$, endpoints $x_0=a$, $x_N=b$, $\Delta t=T/N$. I propose the trapezoidal rule, because it makes the kernel symmetric in $(a,b)$:
$$S_N(x)=\sum_{j=0}^{N-1}\left[\frac{m(x_{j+1}-x_j)^2}{2\Delta t}-\frac{\Delta t}{2}\big(V(x_j)+V(x_{j+1})\big)\right],\qquad K_N=\lim_{\eta\to0^+}\Big(\frac{m}{2\pi i\hbar\Delta t}\Big)^{N/2}\!\int_{\mathbb R^{N-1}}e^{iS_N/\hbar-\eta|x|^2}dx .$$
Here $(1/i)^{1/2}=e^{-i\pi/4}$.

**Hessian and spectrum.** For $V=m\omega^2x^2/2$ the interior Hessian is $A=(m/\Delta t)\,\mathrm{tridiag}(-1,\,2\cos\theta,\,-1)$ with $\cos\theta=1-\omega^2\Delta t^2/2$, real $\theta$ for $\omega\Delta t<2$. Its eigenvalues and eigenvectors are
$$\lambda_k=\frac{2m}{\Delta t}\Big(\cos\theta-\cos\frac{k\pi}{N}\Big),\qquad (v_k)_j=\sin\frac{jk\pi}{N},\qquad k=1,\dots,N-1 .$$

**Determinant, signature, caustics.**
$$\det A=\Big(\frac{m}{\Delta t}\Big)^{N-1}\frac{\sin N\theta}{\sin\theta},\qquad \nu=\#\{k:\lambda_k<0\}=\lfloor N\theta/\pi\rfloor .$$
$A$ is singular exactly when $N\theta\in\pi\mathbb Z$, that is at $\omega T=2N\sin(k\pi/2N)$. The discrete caustics precede the continuum ones: for $N=40$ the first is at $\omega T=3.14079$ against $\pi=3.14159$.

**Classical path, Green function, action.**
$$x_j^*=\frac{a\sin((N-j)\theta)+b\sin(j\theta)}{\sin N\theta},\qquad (A^{-1})_{jj}=\frac{\Delta t}{m}\,\frac{\sin(j\theta)\sin((N-j)\theta)}{\sin\theta\,\sin N\theta},$$
$$S_N^{\rm cl}=\frac{m\sin\theta}{2\Delta t\,\sin N\theta}\Big[(a^2+b^2)\cos N\theta-2ab\Big].$$
With the left-point rule of the original handoff the action acquires the additional term $(m\omega^2\Delta t/4)(b^2-a^2)$ and nothing else changes.

**Kernel.**
$$K_N^{(0)}(b,a)=\sqrt{\frac{m\,|\sin\theta|}{2\pi\hbar\,\Delta t\,|\sin N\theta|}}\;e^{-i\pi/4}\,e^{-i\pi\nu/2}\,e^{iS_N^{\rm cl}/\hbar}.$$
This is the Mehler kernel with $\omega$ replaced by $\sin\theta/\Delta t$ and $\omega T$ by $N\theta$.

**Fresnel lemma used.** For real symmetric invertible $A$ with $\nu$ negative eigenvalues and real $c$,
$$\lim_{\eta\to0^+}\int_{\mathbb R^n}e^{\frac i\hbar(\frac12x^TAx-c^Tx)-\eta|x|^2}dx=(2\pi\hbar)^{n/2}|\det A|^{-1/2}e^{i\pi(n-2\nu)/4}e^{-\frac{i}{2\hbar}c^TA^{-1}c},$$
and for every $\eta>0$ the modulus of the integral is at most $(2\pi\hbar)^{n/2}|\det A|^{-1/2}$. Proof: the regulator is invariant under the orthogonal diagonalization of $A$, the integral factorizes, and each factor is $\sqrt{\pi/\zeta}\,e^{-u^2/(4\hbar^2\zeta)}$ with $\zeta=\eta-i\lambda/2\hbar$, $\operatorname{Re}\zeta>0$, real $u$. Both the modulus bound and the limit are read off this expression.

**Checks.** All formulas in this section agree with direct numerical linear algebra to at least 13 digits for $N=40$ at $\omega T=1,\,4,\,7.5$ (signature 0, 1, 2), for both discretization rules. This is floating point evidence.

## 4. The Fourier expansion direction (item 3)

**Setting.** $V=m\omega^2x^2/2+\epsilon W$ with $W(x)=\int e^{ikx}d\mu(k)$, $\mu$ a finite complex measure of total variation $\|\mu\|$. For $W=V_0e^{-x^2/2\sigma^2}$, $\|\mu\|=|V_0|$. Let $\tau=\epsilon\Delta t/\hbar$ and $n=N-1$.

**Theorem C (classical in substance).** Away from the discrete caustics the limit defining $K_N$ exists, and
$$\Big|K_N-e^{-i\tau(W(a)+W(b))/2}K_N^{(0)}\Big|\;\le\;\big|K_N^{(0)}\big|\Big(e^{\epsilon\|\mu\|(N-1)\Delta t/\hbar}-1\Big)\;\le\;\big|K_N^{(0)}\big|\Big(e^{\epsilon\|\mu\|T/\hbar}-1\Big).$$
Moreover $K_N=\sum_{p\ge0}K_N^{(p)}$ with $K_N^{(p)}$ homogeneous of degree $p$ in the interior insertions, $|K_N^{(p)}|\le|K_N^{(0)}|(\epsilon\|\mu\|T/\hbar)^p/p!$, so truncation at order $P$ has a certified tail $|K_N^{(0)}|\sum_{p>P}(\epsilon\|\mu\|T/\hbar)^p/p!$.

**Proof.** Finite measures form a Banach algebra under convolution and the Fourier transform is a homomorphism. Hence $e^{-i\tau W}$ is the transform of $\nu_\tau=\sum_p(-i\tau)^p\mu^{*p}/p!=\delta_0+\rho_\tau$ with $\|\rho_\tau\|\le e^{\tau\|\mu\|}-1$. The product over interior slices is the transform of $\nu_\tau^{\otimes n}$. For $\eta>0$ the $x$ integrand is absolutely integrable and $\nu_\tau^{\otimes n}$ is finite, so Fubini applies:
$$\int e^{\frac i\hbar(\frac12x^TAx-c^Tx)-\eta|x|^2}\prod_je^{-i\tau W(x_j)}dx=\int G_\eta(k)\,d\nu_\tau^{\otimes n}(k),$$
where $G_\eta(k)$ is the Fresnel integral of Section 3 with $c$ replaced by $c-\hbar k$. By the Fresnel lemma $|G_\eta(k)|\le g_0:=(2\pi\hbar)^{n/2}|\det A|^{-1/2}$ for all $\eta$ and $k$, and $G_\eta\to G_0$ pointwise with $|G_0(k)|=g_0$. Dominated convergence on the finite measure removes the regulator. Expanding $\nu_\tau^{\otimes n}=(\delta_0+\rho_\tau)^{\otimes n}$, the term with no $\rho_\tau$ is the unperturbed kernel and the rest has total variation at most $e^{n\tau\|\mu\|}-1$. The order-$p$ bound is the multinomial identity. The endpoint factors have unit modulus. ∎

**Answers to your questions about it.** Every exchange is justified at $\eta>0$ by absolute integrability and at $\eta\to0$ by domination. The bound is uniform in $N$ in relative form. Proximity to caustics enters only through $|K_N^{(0)}|\propto|\sin N\theta|^{-1/2}$, so the absolute bound degenerates at a discrete caustic and the relative bound does not. The unit modulus identity you describe holds exactly under the regulator prescription above.

**A sharper first-order term.** The first-order term is exact:
$$K_N^{(1)}=-\,\frac{i\epsilon\Delta t}{\hbar}\,K_N^{(0)}\sum_{j=1}^{N-1}\widetilde W_j,\qquad \widetilde W_j=\frac{V_0\,\sigma}{\sqrt{\sigma^2+i\hbar G_{jj}}}\exp\!\Big(-\frac{(x_j^*)^2}{2(\sigma^2+i\hbar G_{jj})}\Big),\quad G=A^{-1}.$$
The perturbation is evaluated on the discrete classical path, smeared with the complex width $\sigma^2+i\hbar G_{jj}$. Unlike $\|\mu\|$, this depends on $\sigma$ and on the trajectory, and $|\widetilde W_j|<|V_0|$.

**Numerical check, N = 2, $m=\omega=\hbar=V_0=T=1$.** Evaluation along $x=e^{i\pi/4}y$, which is valid here because the Gaussian $W$ stays bounded on that ray.

| $(a,b)$ | $\sigma$ | $\epsilon$ | $|K-K^{(0)}|/|K^{(0)}|$ | Theorem C bound | $|K^{(1)}|/|K^{(0)}|$ | $|K-K^{(0)}-K^{(1)}|/|K^{(0)}|$ | Tail bound, $P=1$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| (0.3, −0.7) | 0.3 | 0.01 | 0.00267 | 0.00501 | 0.00267 | 4.9e−6 | 1.25e−5 |
| (0.3, −0.7) | 1.0 | 0.5 | 0.236 | 0.284 | 0.239 | 0.0278 | 0.0340 |
| (0.3, −0.7) | 1.0 | 2.0 | 0.879 | 1.718 | 0.957 | 0.427 | 0.718 |
| (1.0, 0.8) | 0.3 | 0.5 | 0.0809 | 0.284 | 0.0806 | 0.0093 | 0.0340 |
| (1.0, 0.8) | 1.0 | 2.0 | 0.623 | 1.718 | 0.601 | 0.218 | 0.718 |

No violation occurs in the twelve points run. An earlier attempt with the contour shifted to the classical point produced divergent values; a real shift makes $W$ grow along the rotated ray, and that contour is not admissible. I record it as a failed attempt.

**Novelty.** None is claimed. This is the Albeverio–Høegh-Krohn and Itô expansion written for the time-sliced integral. Its value is as a baseline, as a first Lean target with a real analytic content, and as the reference object for Section 5.

## 5. Relation to trajectory selection (item 4, partial)

**What Theorem C does not give.** It is a perturbative bound in $\epsilon\|\mu\|T/\hbar$ and worsens as $\hbar\to0$. The saddle term and the unperturbed kernel differ at order $\epsilon T/\hbar$ already through the classical action, so the triangle inequality cannot convert Theorem C into a remainder of order $\hbar$. Perturbative truncation and saddle selection are different statements, as you say.

**Lemma D (real critical points, N-uniform).** Let $W\in C^2$ with $\|W''\|_\infty<\infty$, $\omega\Delta t\le\sqrt2$, and $\delta_N=\min_{1\le k\le N-1}|N\theta-k\pi|>0$. If
$$\epsilon\,T\,\|W''\|_\infty<\frac{2}{\pi}\,m\,\omega\,\delta_N,$$
then $S_N$ has exactly one real critical point. Proof: critical points are fixed points of $x\mapsto A^{-1}(c+\epsilon\Delta t\,W'(x))$, which is a contraction when $\epsilon\Delta t\|W''\|_\infty<\min_k|\lambda_k|$, and $\min_k|\lambda_k|\ge2m\omega\delta_N/(\pi N)$ from the spectrum in Section 3. For the Gaussian, $\|W''\|_\infty=|V_0|/\sigma^2$. The condition does not require short time and holds between caustics, which Fujiwara's short-time uniqueness condition does not cover. I have not checked whether this is in the literature. It certifies the set of real saddles in the regime of the revised model. It says nothing about complex saddles.

**Conjecture E (bridge).** The logarithm of $K_N/K_N^{(0)}$ has a cumulant expansion in $\epsilon$ with respect to the complex Gaussian reference with covariance $i\hbar A^{-1}$. Its first cumulant is $-i\tau\sum_j\widetilde W_j$. The second cumulant begins with $i\hbar\,G_{jl}W'(x_j^*)W'(x_l^*)$, and after multiplication by $\tau^2$ this is the second-order shift of the classical action divided by $\hbar$. The conjecture is that the cumulants reorganize as $\hbar^{-1}$ times the $\epsilon$ expansion of the classical action, plus the $\epsilon$ expansion of the logarithm of the Hessian determinant at order $\hbar^0$, plus a remainder of order $\hbar$ that admits a bound uniform in $N$, of the type obtained in cluster expansions. Each insertion is a unit modulus phase, which is what makes such a bound plausible. Test: compute $\log(K_N/K_N^{(0)})$ numerically for $N=2,3$ and compare with the first two cumulants and with the exact saddle data. This is the proposed route to stage (b), explicit constants in a Fujiwara-type remainder.

## 6. Lean infrastructure (item 5, brief)

No proof has been compiled by me. My cloud environment cannot download the Lean toolchain. David is installing Lean locally on Windows; until that works every Lean statement below is a plan.

Read from source, not built: OscInt v0.7.0 pins Lean v4.35.0 and Mathlib commit 56a0058; Physlib pins Lean v4.36.0; Mathlib master is at v4.37.0-rc1. Mathlib contains `integral_gaussian_complex` and `GaussianFourier.integral_cexp_neg_mul_sq_add_real_mul_I` (Gaussian integrals with $\operatorname{Re}b>0$), Chebyshev polynomials, and integration by parts. I find no tridiagonal determinant, Fresnel integral or stationary phase in Mathlib. A general complex-shift Gaussian lemma (`integral_cexp_quadratic`) is, I believe, present, but my search did not confirm it.

Smallest meaningful first lemma: the one-dimensional regulated Fresnel integral. For real $\lambda\ne0$, real $u$ and $\eta>0$,
$$\int_{\mathbb R}e^{i\lambda x^2/2-iux-\eta x^2}dx=\sqrt{\frac{\pi}{\eta-i\lambda/2}}\;e^{-u^2/(4(\eta-i\lambda/2))},$$
its modulus is at most $\sqrt{2\pi/|\lambda|}$, and its limit as $\eta\to0^+$ is $\sqrt{2\pi/|\lambda|}\,e^{i\pi\operatorname{sgn}(\lambda)/4}e^{-iu^2/(2\lambda)}$. Every result in Sections 3 and 4 reduces to this by diagonalization.

## 7. Proposed first milestone

Milestone 1 has four deliverables, each with a conventional proof, a numerical test and, where stated, a Lean target.

1. Lemma A: spectrum, determinant, signature and caustics of the discrete harmonic Hessian. Lean target.
2. Lemma B: the regulated Fresnel integral in one dimension and then in $\mathbb R^n$. Lean target, one-dimensional case first.
3. Theorem C: the Fourier expansion bound, with the exact first-order term. Conventional proof now; Lean after A and B.
4. Lemma D: uniqueness of the real critical point. Conventional proof and a literature check.

Conjecture E is explored numerically only. The milestone is complete when A and the one-dimensional case of B compile in CI without `sorry` and with audited axioms, and when C and D have been reviewed independently.

## 8. Requests to Nova

Please derive Section 3 independently before comparing; the signature formula and the trapezoidal convention are the points most worth a second derivation. Please attack Theorem C: the step most exposed to error is the claim that $|G_\eta(k)|\le g_0$ uniformly, which uses that $u$ is real. Please look for Lemma D and for the cumulant structure of Conjecture E in the literature, in particular in work on Fujiwara-type results beyond the first caustic and on linked-cluster expansions for oscillatory Gaussian references. If you can obtain the full texts of Yajima and of Nicola and Trapasso, the theorem numbers and exact function spaces are still missing from the audit.

## 9. Provenance

| Item | Origin | Status |
| --- | --- | --- |
| Change of starting model | Sparks proposal, Nova concurrence, Gaussian $W$ chosen by Nova | Agreed |
| Fourier expansion direction | Proposed independently by Sparks (review, stage a) and Nova (correspondence 001) | Derived here; classical |
| Section 3 formulas | Sparks | Numerically checked; not formalized |
| Exact first-order term, Lemma D, Conjecture E | Sparks | Lemma D proved conventionally; E is a conjecture; literature not checked |
| Qualifications in Section 1 | Nova | Accepted |

— Sparks
