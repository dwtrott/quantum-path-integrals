# Research correspondence 006: Sparks → Nova and David
**Date:** November 23, 2026
**Project:** Certified Trajectory Selection in Quantum Path Integrals
**Subject:** Novelty audit of T. Hargé's papers against the connected-tree estimate (correspondence 003, Section 5) and Proposition G (correspondence 004, Section 4)
**Status:** Literature audit from the arXiv TeX sources supplied by David. No new mathematics.

## 0. Verdict

The method is not new. Hargé's paper arXiv:1305.3738 contains, for the Euclidean continuum heat kernel, the same chain of steps that Nova proposes in correspondence 003 and that I extend in correspondence 004: a Dyson expansion written through the Fourier measure of the perturbation, a formal logarithm, the Brydges–Kennedy tree graph equality to restore convergence, a smallness condition on an exponential Fourier moment, and the conclusion that the logarithm of the kernel is a classical term divided by $h$ plus a correction regular in $h$. His worked example is the Gaussian perturbation $c(x)=ke^{-x^2/2}$ with $|k|$ small, which is our model.

What his papers do not contain, in the parts I read, is the real-time statement in $\hbar$ across caustics, the time-sliced statement uniform in $N$, and explicit constants. Those are the only places where a contribution could remain, and each needs a further check, in particular against his paper [Ha7], which I do not have.

The connected-tree estimate and Proposition G should therefore be described as a transposition of Hargé's method to the real-time, time-sliced setting, and not as a new technique.

## 1. What I read

| Paper | Read | Not read |
| --- | --- | --- |
| T. Hargé, *Borel summation of the semi-classical expansion of the partition function associated to the Schrödinger equation*, arXiv:1305.3738v1 (May 2013) | Introduction; Section 2 (notation and main results); Section 3.1 (tree graph equality and the remarks after it); Section 4 through the proposition defining $\phi$ and $w$ | Borel summation sections, Gaussian summation, Morse lemma, proofs, appendices |
| T. Hargé, *Borel summation of the small time expansion of the heat kernel. The scalar potential case*, arXiv:1301.7742v1 (January 2013), cited there as [Ha4] | Part of the introduction; the main theorems of Section 3 | Proofs |
| T. Hargé, *A deformation formula for the heat kernel*, arXiv:1302.1810v1, cited as [Ha6] | Introduction; the main theorem | Proofs |

The sources carry labels and not printed numbers, and I did not compile them. I therefore cite results by section and description. Printed theorem numbers remain to be filled in from the PDFs.

## 2. Correspondence between Hargé's objects and ours

Hargé studies $H=-h^2\partial_x^2+\tfrac14(\omega x)^2-c(x)$ on $\mathbb R^\nu$ with $c(x)=\int e^{ix\cdot\xi}d\mu(\xi)$, and the kernel of $e^{-tH/h}$ for real $t\in(0,T)$ and $\operatorname{Re}h>0$.

| Hargé (1305.3738) | This project |
| --- | --- |
| Deformation formula, Section 4: $v=1+\sum_mv_m$, with $v_m=(t/h)^m\int\!\!\int e^{-h\,\Omega_t.\xi\otimes\xi}e^{iq_t^m(s)\cdot\xi}\,d\mu^{\otimes}d^ms$, where $q_t$ is the harmonic classical path and $\Omega_t$ is built from the harmonic Dirichlet Green function | Fourier expansion of correspondence 002, Section 4: $Z=\int G_0(k)/G_0(0)\,d\nu^{\otimes n}$ with $e^{ik\cdot x^*-i\hbar k^TGk/2}$ |
| Remark after the deformation formula: $|v_m|\le\frac1{m!}(t/h)^m\|\mu\|^m$, using $\Omega_t\ge0$ | Theorem C: $|Z^{(p)}|\le(\epsilon\|\mu\|T/\hbar)^p/p!$, using unit modulus |
| Tree graph equality, Section 3.1, attributed to Brydges–Kennedy (J. Stat. Phys. 48, 1987) and Brydges' Lausanne notes: $A=e^B$ with $B$ a sum over trees with interpolation parameter $\theta_{j,k,g}$ equal to the minimum along the tree path | Nova's forest interpolation, correspondence 003, Section 5.2, with the same minimum-along-path weights |
| Section 3.1, following proposition: a mould identity of the form $C_h=\exp(\tfrac1hE+R_h)$ | The splitting $\log Z=\alpha/\hbar+\beta+R_N$ of Proposition G |
| Section 4: $Q_m$ defined as tree sums without the interpolation exponent; $\varphi=\sum t^{2m-2}Q_m$; kernel $=(\text{harmonic prefactor})\times e^{-\phi/h+w}$ under $4T^2M_\mu<1$, $M_\mu=\int e^{5|\xi|}d|\mu|$ | Candidate theorem of correspondence 004: $K_N=K_N^{(0)}\cdot|\det A/\det A_\epsilon|^{1/2}e^{i(S_\epsilon-S_0)/\hbar}(1+r_N)$ under $(27/4)|\epsilon|TL_{r+s}M/r^2<1$ |
| Remark in Section 4: $\phi$ solves the Hamilton–Jacobi equation $\partial_t\phi+(\partial_x\phi)^2=\tfrac14(\omega x)^2-c(x)$ | Identification of $\alpha$ with the classical action |
| Section 2, proposition on the heat kernel: under $4T^2e^TM_{\mu,\varepsilon}<1$ the kernel equals $h^{-\nu/2}e^{-\phi/h}$ times a Laplace transform in $h$, that is, the $h$-expansion of $e^w$ is Borel summable | Stronger than our order-$\hbar$ remainder: all orders, with resummation |
| Section 2, theorem: Borel summability of the partition function, with further assumptions | Not an object of this project |

Two details confirm how close the settings are. His smallness condition has the form $T^2\times$(exponential Fourier moment)$<$ constant; ours has the same form in the free case, where $M\le T/(4m)$. His introduction states that the potential with $c(x)=ke^{-x^2/2}$ and small $|k|$ has no real critical point other than the origin and that tunnelling is not considered, which is the single-real-saddle regime of Lemma D.

## 3. What the papers say about real time

In 1305.3738 the time $t$ is real and positive and $h$ lies in the right half plane. A remark in Section 4 states that treating complex $h$ and complex $t$ together certainly needs an extra assumption, and points to [Ha7, Prop. 4.5, case 2]. Convergence of the deformation formula for arbitrary $t$ there uses the positivity of $\Omega_t$. Real time corresponds to imaginary $t$, so the semiclassical statement of that paper does not cover the Schrödinger propagator as written.

In 1301.7742 the potential is $\alpha x^2-c(x)$ with $\omega\in\mathbb R\cup i\mathbb R$, and the introduction states that the method treats complex time with $\operatorname{Re}t\ge0$. The main theorem for the harmonic case gives, for small $|t|<T$ in the half disk, an expansion in powers of $t$ with remainders bounded by $K\,r!\,\kappa^{-r}|t|^r$, uniformly as $\operatorname{Re}t\to0^+$, and the following remark says that this yields information about the Schrödinger kernel. This is a small-time expansion, not a semiclassical one, and it is restricted to small times.

In 1302.1810 the main theorem gives the deformation formula for a time-dependent setting, with the conjugate kernel smooth up to imaginary times in an interval $i(-T_c,T_c)$. Again the time is small.

So the real-time deformation formula, which is the continuum counterpart of Theorem C, is established for small times in these papers, in addition to the older work of Itô and of Albeverio and Høegh-Krohn that Hargé cites for it.

## 4. Where a contribution could remain

Each item below is a possibility, not a claim.

1. **Semiclassical remainder in real time.** A bound of order $\hbar$ on the remainder of the real-time kernel after the action and determinant terms. Hargé's semiclassical paper stops at real $t$; whether [Ha7] or later work closes this is unknown to me.
2. **Validity between and beyond caustics.** Our estimates hold for any $T$ with $|\sin N\theta|\ge s_0$ and carry the Maslov phase. Hargé's real-time statements are for small times. In real time the unit modulus of the oscillatory factor replaces the positivity he uses, which is the one structural difference between the two settings.
3. **Time-sliced statement uniform in $N$.** Hargé works in the continuum. A finite-$N$ theorem with constants independent of $N$ is a different object, relevant to Feynman's definition and to the Fujiwara line of work.
4. **Explicit constants.** His constants $K$, $\kappa$, $\alpha_\nu$ are existential. Ours are computable.
5. **Formal verification.** Not a mathematical novelty.

My assessment is that item 1 alone would probably be regarded as accessible from Hargé's method by a specialist. Items 2, 3 and 4 together could support a modest paper that is explicitly positioned as an extension of Hargé to real time and to the time-sliced integral. None of this reaches the original objective, certified selection among several saddles, which Hargé's hypotheses also exclude.

## 5. What we can take from Hargé

- His route to the identification of the classical term is different from mine: he shows that $(t/h)\varphi+w$ satisfies a nonlinear equation, so that its exponential solves the kernel equation, and reads off the Hamilton–Jacobi equation for $\phi$. In the continuum this is cleaner than matching asymptotics. At finite $N$ the uniqueness argument of correspondence 004 still seems the simpler one.
- His mould proposition $C_h=\exp(E/h+R_h)$ is the algebraic form of the splitting in Proposition G and is a ready-made check on it.
- His choice of norms, with weight $e^{5|\xi|}$ and for Borel summability $e^{\varepsilon\xi^2+5|\xi|}$, indicates which moments are needed at each stage.
- The references [Fu-Os-Wi] (Fujiwara, Osborn, Wilk, Phys. Rev. A 25 (1982) 14–34) and [On] (Onofri, Am. J. Phys. 46 (1978) 379–382) are cited for the heuristic form of the same tree expansion and belong in our literature list.

## 6. Open items

1. Obtain [Ha7]: T. Hargé, *Some remarks on the complex heat kernel on $\mathbb C^\nu$ in the scalar potential case* (2013). Its Proposition 4.5 is the place where complex $t$ and $h$ are treated together.
2. Search for later papers citing Hargé 1305.3738 that treat the real-time kernel or oscillatory integrals. One known descendant is Albeverio, Cangiotti and Mazzucchi, arXiv:1810.13158, on small-time expansions for stochastic differential equations.
3. Fill in the printed theorem numbers from the PDFs.
4. Decide the direction (for David): pursue the extension described in Section 4 as a bounded first paper, or move directly to the first two-saddle configuration, where neither Hargé's hypotheses nor our tree domain apply.

## 7. Recommended change to the record

The ledger entries for the connected-tree bound and Proposition G should carry the note: "method due to Hargé (arXiv:1305.3738) in the Euclidean continuum setting; possible novelty limited to real time, caustics, finite $N$ and explicit constants; [Ha7] unchecked."

— Sparks
