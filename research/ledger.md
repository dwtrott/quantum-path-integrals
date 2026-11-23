# Contribution ledger

| Date | Item | Contributor | Reviewer | Status |
| --- | --- | --- | --- | --- |
| 2026-11-23 | Change of starting model from the quartic oscillator to the harmonic oscillator with a bounded perturbation | Sparks (proposal), Nova (concurrence, Gaussian perturbation) | David | Agreed in correspondence |
| 2026-11-23 | Discrete harmonic kernel: spectrum, determinant, signature, caustics | Sparks | Nova (pending) | Conventional derivation; numerically checked; not formalized |
| 2026-11-23 | Fourier expansion bound (Theorem C), classical in substance | Sparks and Nova, independently proposed | pending | Conventional proof; checked numerically at N = 2; not formalized |
| 2026-11-23 | Uniqueness of the real critical point (Lemma D) | Sparks | pending | Conventional proof; literature not checked |
| 2026-11-23 | Cumulant bridge (Conjecture E) | Sparks | pending | Conjecture |
| 2026-11-23 | Connected-tree cumulant bound, exponential-moment extension, zero-free region | Nova | Sparks | Reconstructed independently; numerically tested, no violation. Method due to Hargé (arXiv:1305.3738) in the Euclidean continuum setting; possible novelty limited to real time, caustics, finite N and explicit constants; [Ha7] unchecked |
| 2026-11-23 | Green function formula, entrywise bound, max-norm uniqueness criterion | Nova | Sparks | Numerically confirmed |
| 2026-11-23 | Order-hbar remainder (Proposition G); identification by uniqueness of asymptotics | Sparks, on Nova's tree representation | Nova (pending) | Derived; unreviewed. Same novelty note as the tree bound (see correspondence 006) |
| 2026-11-23 | Conjecture F (zero-free radius and saddle degeneracy) | Sparks | pending | Numerical observation |
| 2026-11-23 | Lean toolchain smoke test | Sparks | CI run 45263613230 | Compiled; standard axioms only |
| 2026-11-23 | Multi-saddle proposal: cusp structure, selection of one complex saddle, reduction to a soft mode | Sparks | Nova (pending) | Proposal with numerical evidence; prior art unchecked |
| 2026-11-23 | Novelty audit of Hargé 1305.3738, 1301.7742, 1302.1810 | Sparks | Nova (pending) | Done from TeX sources, partial reading; see correspondence 006 |
| 2026-11-23 | Lean: closed form of the discrete harmonic determinant recurrence (`scaledDet_mul_sin`, `scaledDet_eq`, `scaledDet_eq_zero_iff`) | Sparks | CI run 45268235546 | Compiled on Lean 4.35.0; standard axioms only; link to `Matrix.det` not yet proved |
| 2026-11-23 | Ha7 follow-up to the Hargé audit (arXiv:1302.3512, Prop. 4.5 case 2) | Nova | Sparks (pending) | Reported in `research/NOVA_RESPONSE_TO_HARGE_AUDIT.md`; novelty of real-time finite-N explicit constants unresolved |
| 2026-11-23 | Proposition H: exact classification of three real stationary points in a symmetric Gaussian well | Nova | Sparks | Conventional proof; each step checked by Sparks without finding an error; numerically confirmed for N = 2, 3, 4, 6, 9 (`prop_h_check.py`); novelty unchecked; not in Lean |
| 2026-11-23 | Selection certificate in one variable (Lemma 1, Lemma L, Theorem S) with a worked example | Sparks | Nova (pending) | Conventional proofs, unreviewed; example hypotheses evaluated in floating point, not certified; expected to be known in substance |
| 2026-11-23 | Review of contour lemma, local Gaussian error constants, polygonal alternative and three-real-saddle contour existence (009) | Nova; original certificate Sparks | Sparks (pending) | Conventional review and numerical evidence; no interval or Lean certification; novelty unchecked |
