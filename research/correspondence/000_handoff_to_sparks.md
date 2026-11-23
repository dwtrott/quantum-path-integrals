# Research Handoff to Sparks (Claude)
## Certified Trajectory Selection in Quantum Path Integrals
**Date:** November 23, 2026  
**Collaborators:** David (human research lead), Nova (ChatGPT), Sparks (Claude)  
**Status:** Project initiation; no novel theorem or verified project proof yet.

## 1. Purpose and shared vision

Sparks, David and Nova have been discussing a collaborative mathematical physics research project. The goal is not to solve a Millennium Prize problem or merely reproduce textbook results. We want to discover a tractable, genuinely new mathematical result concerning Feynman's **real-time path integral formulation of quantum mechanics**, using AI-assisted conjecture generation, numerical experimentation and Lean 4 formal verification. We want a result that could eventually withstand peer review and potentially be published.

The motivating question is whether one can identify the classical and complex saddle contributions that matter to a quantum transition, with mathematically justified bounds on the contributions omitted. The aim is not to assert that particles literally traverse selected trajectories or to settle competing interpretations of quantum mechanics. Mathematical explainability, explicit assumptions and predictive accuracy are the goals.

This is a new research direction in an ongoing human–AI collaboration, not our first research project together. David has a background in physics and applied mathematics. He intends to concentrate on orchestration, independent scrutiny, research direction and translating formal results into understandable mathematics. Sparks and Nova should make substantive independent intellectual contributions and challenge one another rather than simply implement David's ideas.

## 2. Mathematical starting point

For fixed endpoints $x(0)=a$ and $x(T)=b$, the formal propagator is

$$K(b,T;a,0)=\int_{a}^{b}\mathcal D x\,\exp(iS[x]/\hbar),\qquad S[x]=\int_0^T [m\dot x^2/2-V(x)]\,dt.$$

The real-time integral is oscillatory, not an ordinary probability measure on paths. Stationarity of the action gives classical trajectories, but real stationary paths alone do not in general capture tunneling or other nonperturbative effects. Complex saddles, contour geometry, interference and Stokes phenomena may be essential.

Start with the one-dimensional quartic anharmonic oscillator, $V(x)=m\omega^2x^2/2+\lambda x^4$, with $m>0$, $\omega>0$, $\lambda\geq0$. Use the harmonic oscillator ($\lambda=0$) as a baseline, avoiding caustic times initially. The double well and Rosen–Morse barrier are follow-on benchmarks.

For $N$ time intervals of length $\Delta t=T/N$, let $x_0=a$, $x_N=b$ and define a consistent time-sliced action, for example

$$S_N(x)=\sum_{j=0}^{N-1}\left[\frac{m(x_{j+1}-x_j)^2}{2\Delta t}-\Delta t\,V(x_j)\right].$$

The corresponding discretized propagator has normalization $(m/(2\pi i\hbar\Delta t))^{N/2}$ times the oscillatory integral over $\mathbb R^{N-1}$. Specify a regulator, e.g. $e^{-\eta\|x\|^2}$ for $\eta>0$, and distinguish fixed-regulator results from removal of the regulator and the continuum limit. Complex contour deformation requires separate analyticity, convergence and boundary justification. Do not silently exchange limits.

## 3. Candidate research directions (not established novelty)

### A. Certified nonstationary contribution bounds

For a smooth compactly supported amplitude $a$ supported where $|\nabla S|\geq\delta>0$, use integration by parts with $v=\nabla S/|\nabla S|^2$:

$$\left|\int a(x)e^{iS(x)/\hbar}dx\right|\leq\hbar\int |\nabla\cdot(a v)|dx.$$

This is a standard starting estimate, not a new discovery. Repeated integration by parts is possible under suitable regularity. Investigate explicit constants for the discretized action, dependence on $N$, and adaptive partitions that isolate near-stationary regions. A bound that becomes useless as $N$ grows must be reported as such.

### B. Certified saddle truncation

If a justified contour decomposition gives $I=\sum_j n_j I_j$, and certified bounds $|I_j|\leq B_j$ are available, then for a retained subset $J$,

$$|I-\sum_{j\in J}n_j I_j|\leq\sum_{j\notin J}|n_j|B_j.$$

For approximate retained contributions add their individual error bounds; for infinite sums establish convergence and a tail estimate. This triangle inequality is elementary; the substantive challenge is obtaining useful certified $n_j$, $B_j$, and error estimates without already solving the full integral. Explore bounds from thimble geometry, validated numerics and resurgence. Absolute errors are primary because destructive interference can make relative error ill-conditioned. Treat Stokes transitions and degenerate saddles explicitly; do not select saddles by action magnitude alone.

### C. Cross-disciplinary discovery

Systematically search for transferable techniques from asymptotic analysis, Morse theory, algebraic topology, complex analysis, validated numerical analysis, dynamical systems and optimization. The motivating analogy is the use of the Euler beta function in Veneziano's amplitude: mathematical structures can acquire new uses long after their discovery. For each candidate analogy, specify an actual mathematical mapping, assumptions, possible failure cases and a small test. Avoid superficial analogies and claims of novelty without literature checks.

## 4. Literature anchors

Bibliographic details and abstracts for these two recent papers were checked on November 23, 2026; full technical results still require independent reading and reproduction.

1. Job Feldbrugge, *Real-Time Path Integral Resurgence: Which Classical Path Contributes?*, preprint, November 9, 2026. https://arxiv.org/abs/2611.12510 — uses high-order saddle expansions and Borel–Padé analysis to infer saddle relevance, with a Rosen–Morse barrier benchmark. Prioritize identifying assumptions, numerical versus rigorous claims, and any existing error estimates.
2. Aurélien Dersy and Matthew D. Schwartz, *Beyond the Dilute Instanton Gas: Resurgence with Exact Saddles in the Double Well*, Phys. Rev. Lett. 137, 181601, October 27, 2026. https://doi.org/10.1103/hg9t-ypvf ; https://arxiv.org/abs/2604.14279 — finite Euclidean-time exact saddles and resurgence for nonperturbative energy splittings. Distinguish its Euclidean setting from our real-time objective.
3. Expand backward through their references on Picard–Lefschetz theory, stationary/nonstationary phase, complex semiclassical trajectories, resurgence and rigorous path integral constructions. Search forward for work that already establishes proposed bounds.

Do not treat any claimed recent AI mathematical breakthrough as a premise for this project. The research must stand on its own.

## 5. Lean and computational infrastructure

- Lean 4 and Mathlib: https://github.com/leanprover-community/mathlib4
- Physlib (physics formalization): https://github.com/leanprover-community/physlib
- A project referred to in discussion as LeanOscillatory was suggested as a possible source of oscillatory integral lemmas. Its exact repository, coverage and compatibility have not been independently verified in this handoff; locate and verify it before depending on it. Search Mathlib and other packages for existing integration-by-parts and stationary phase results.
- Begin with the smallest useful formal theorem (e.g. a one-dimensional compact-support nonstationary phase bound), then generalize to higher dimensions and the concrete action. Do not attempt to formalize all quantum mechanics first.
- Install Lean via elan and use Lake. Pin a mutually compatible Lean version and dependency commits; commit `lean-toolchain`, Lake configuration and `lake-manifest.json`. Use cached Mathlib builds where possible. Run `lake build` locally and in GitHub Actions.
- An actual compiler feedback loop is essential: inspect goals and errors, revise, and check again. Local execution on David's machine is preferred for fast iteration when access is authorized; GitHub CI is independent reproducibility verification. Both AI systems' access must be established, not assumed.
- No final theorem may depend on `sorry`, `admit`, or unapproved axioms. Inspect `#print axioms` on named results, check theorem statements against the mathematical claim, and audit trusted dependencies. A successful compilation alone is insufficient.
- Use Python (NumPy/SciPy, mpmath, SymPy as appropriate) for exploratory numerics, exact/spectral reference calculations, and high-precision comparisons. Use interval arithmetic or validated bounds when claiming numerical certification. Ordinary floating-point agreement is evidence, not proof.

## 6. Shared GitHub workflow

Proposed private repository: `quantum-path-integrals` (not yet created or configured as part of this handoff).

```text
README.md
AGENTS.md                    # shared research rules for AI collaborators
research/
  literature.md
  journal.md
  conjectures/
  counterexamples/
proofs/                      # conventional, human-readable arguments
QuantumPathIntegrals/         # Lean source modules
experiments/                  # scripts, configs, reproducible outputs
reports/                      # benchmark and proof-review reports
.github/workflows/            # Lean build and numerical tests
lean-toolchain
lakefile.toml                 # or lakefile.lean
lake-manifest.json
```

Use issues for scoped questions, separate branches/worktrees for Nova and Sparks, pull requests for review, and required CI checks before merging accepted proofs. Suggested labels: `literature`, `conjecture`, `counterexample`, `formalization`, `experiment`, `verified`, `blocked`. Record the environment, parameters, random seeds if relevant, commands and commit hashes with each result. Keep failed attempts, but clearly label them. Do not include credentials or unrelated private material.

David coordinates task assignment and acceptance. Nova and Sparks should periodically work independently on identical questions before seeing one another's solution. Neither should impersonate the other or claim a compiler run that did not occur. No assumption of persistent autonomous execution: use explicitly established agent environments and task triggers.

## 7. First milestones and acceptance criteria

1. **Literature and infrastructure audit.** Read the two anchors; produce a table of existing results, hypotheses, limitations and candidate gaps, with page/section references. Inventory actual Lean lemmas and compatible packages. Confirm local and CI access with a small proof and record the build result.
2. **Baseline reproduction.** Implement the discretized harmonic oscillator and compare with the exact propagator away from caustics, using a consistent phase convention. Document regulator, discretization error and numerical error separately. Implement a small quartic case and compare against an independent spectral or numerical reference.
3. **First formal result.** State and prove a standard compact-support nonstationary phase estimate in conventional mathematics and Lean. Include a verified example and a failure case when the phase derivative vanishes.
4. **Candidate novelty.** Propose two or three narrow extensions, with exact hypotheses, proof obligations, test cases, existing literature overlap, estimated effort and a falsification strategy. Rank by mathematical value and tractability. Do not claim novelty merely because a formalization is absent.
5. **Research result.** Pursue the most promising candidate; produce a readable proof, formal verification of its precise scope, and a reproducible numerical demonstration. State explicitly whether the result concerns fixed dimension, regulator removal or the continuum propagator. Publication is a goal, not a guaranteed outcome.

Use a fixed, agreed token/compute budget per investigation; numerical falsification and library search should precede expensive formal proof attempts. Stop or redirect when assumptions fail, bounds are vacuous, or the literature establishes the result already.

## 8. Scientific attribution and human readability

Every substantial result should have a conventional mathematical statement and proof, a corresponding Lean theorem when feasible, and a concise explanation of mathematical significance and implications. Physical interpretation may remain open. Formal verification establishes a statement under its assumptions, not the truth of a physical model or the novelty of a result.

Keep a contribution ledger with date, issue/commit, contributor (David, Nova, Sparks), proposed insight, proof or counterexample, reviewer and verification status. Credit actual intellectual contributions, including AI-generated ones, transparently. Human authorship and accountability must comply with the eventual venue's policies; do not imply that AI systems qualify as authors under all publication rules.

## 9. Immediate request to Sparks

Please review this handoff critically rather than accept it unconditionally. Return:

- Your understanding of the objective and any mathematical objections or missing assumptions.
- A concise assessment of the two recent papers and how our proposal differs from their established results (or mark what needs full-text review).
- A concrete Lean/local setup and CI plan, including what can actually be run in your environment and how Nova can reproduce it.
- Two or three precise candidate theorem statements, including one deliberately modest starting theorem, expected difficulty, and how to disprove or test each.
- A proposed first computational experiment with independent reference, error metrics and parameter ranges.
- Cross-disciplinary leads that could yield a genuinely useful mathematical transfer.

Do not spend substantial resources attempting a complete path integral formalization yet. First agree on a small, testable theorem and the shared environment. David will coordinate the next round of independent investigation and review.
