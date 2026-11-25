# 052. Sparks: novelty audit of the candidate results

Date: 2026-11-24. Author: Sparks (Claude). Status: literature search. Based on web searches and on
machine-extracted summaries of full texts; I have not read the papers line by line. Absence of a result in
these searches is not proof that it is absent from the literature. Purpose: David's decision on whether the
project continues.

## Results that are known and are not to be claimed

1. The selection rule as such: Picard–Lefschetz theory (023).
2. Propositions F and T: sufficient conditions within the saddle point method. Pemantle and Wilson,
   arXiv:0903.3585, prove expansions for complex-phase integrals in d variables by contour shifting in complex
   d-space, under the hypothesis that the real part of their phase is nonnegative on the chain and vanishes at
   the stationary points. Our propositions are that hypothesis applied on a shifted or tilted plane.
3. Nova's concavity theorem (020). For a function holomorphic on a tube and bounded in the real directions,
   the logarithm of its supremum over real translations is a plurisubharmonic function that depends only on
   the imaginary part, hence convex. The function m(y) of 020 is minus ħ times that logarithm. This is
   classical several complex variables, the multivariable form of the three-lines theorem.
4. The order-ħ coefficient a_1 (048): the standard two-loop term.
5. The type of theorem in Nova 050. Gradinaru and Hagedorn, Numer. Math. 2014, prove convergence of a
   wavepacket time-splitting that separates the local quadratic part of the potential from the remainder,
   with L² error linear in the semiclassical parameter and quadratic in the time step. Their constants are
   not explicit, nothing is certified, and they do not treat overlaps or the parity cancellation.

## A new and close reference

arXiv:2606.28271, "Which Saddles Contribute? The South-East Rule for Multidimensional Integrals" (June 2026).
It gives an algorithm for the relevant complex critical points of integrals with a real analytic phase, by
rotating the phase and tracking Stokes transitions in the Borel plane. According to the extracted summary:
it is not stated as a theorem with proof; adjacency is found numerically by Padé analysis; complex-valued
phases are mentioned only as a case where the method can be inconclusive; path integrals are deferred to
future work; conjugate points, certified arithmetic and error bounds are absent. This is the same question
as ours, asked in 2026, and it leaves open the setting we work in.

## Candidates for which the searches found no prior statement

A. Certified selection. Interval-arithmetic proofs that a named complex saddle is selected: 018 (two
   variables, forbidden target), Nova 034 and 036 (N slices, index zero). The closest works (Shoji and
   Trailović 2510.06334; the South-East rule; numerical steepest descent) are numerical.
B. The continuum criterion. Selection for all large N from the sign of the lowest eigenvalue of a
   Sturm–Liouville problem built from the complex classical path (033, Nova 034 and 035). The coherent state
   propagator literature computes the second variation for the prefactor; I found no use of it as a
   selection criterion.
C. The one-mode surface for a path of Morse index one (038, Nova 039 to 046). Unfinished: the middle region
   is open.
D. Explicit certified constants, uniform in the Strang mesh, for the overlap bound of Nova 050, with the
   parity cancellation that gives order ħ. Modest: constants for a known type of theorem.

Searches not done: exact WKB (Aoki, Kawai, Takei; Delabaere, Howls), Fedoryuk's multidimensional saddle point
method, and the Russian literature on the complex germ (Maslov), which is the most likely home of B.

## Assessment

Nothing here is a new rule for which paths contribute. What appears defensible, if A and B survive the
remaining searches, is a paper on certified selection of complex classical paths in time-sliced wavepacket
path integrals: finite-N certificates, a continuum eigenvalue criterion with proof, and explicit constants.
It is a methods and rigour contribution in a question that the 2026 literature treats numerically.
