# 023. Sparks: prior art on complex saddle selection, first pass

Date: 2026-11-23. Author: Sparks (Claude). Status: literature search only. I have read abstracts and search
excerpts, not the papers. Every characterisation below needs confirmation against the full text before it is
relied on.

## Findings

The question "which complex classical trajectories contribute" is an established research topic with a
standard answer and several practical criteria.

1. Picard–Lefschetz theory gives the necessary and sufficient condition: a saddle contributes exactly when its
   steepest ascent manifold meets the original integration cycle, with an integer intersection number.
   Tanizaki and Koike, arXiv:1406.2386, apply this to the real-time path integral. Feldbrugge and Turok,
   arXiv:2207.12798, treat existence of real-time path integrals in this framework. Nishimura and
   collaborators, arXiv:2307.11199, compute thimbles for discretised real-time tunnelling.
2. Shoji and Trailović, arXiv:2510.06334, give a numerical method for intersection numbers in discretised
   real-time path integrals with up to about 20 variables. This is the closest work to ours in aim.
3. Coherent-state and wavepacket propagation: Tomsovic and collaborators, arXiv:1804.10511, give a search
   technique for the contributing complex saddle trajectories in many degrees of freedom. Wang and Tomsovic,
   arXiv:2107.08799, distinguish exposed and hidden saddles and state a simple criterion based on the real
   reference trajectory. Aguiar and collaborators, Phys. Rev. E 69, 066204, and arXiv:0704.2510, discuss
   noncontributing trajectories and phase space caustics for the coherent state propagator; the latter states
   that a general inclusion criterion exists and is hard to verify in practice.
4. Kicked maps: Shudo and Ikeda with collaborators (arXiv:nlin/0205048, nlin/0301004, nlin/0105067 and the
   review arXiv:2604.12926) study the complex trajectories of the Hénon and standard maps, their Stokes
   geometry, and the removal of trajectories whose imaginary action diverges.

## Consequences for our claims

The selection rule as a statement is not new: it is the intersection number rule. The observations in 019
that smallest Im Φ is not the criterion, and that near-real saddles are the ones that contribute, match what
this literature reports. Propositions F and T are sufficient conditions of the saddle point method and should
be presented as such.

What I have not found in this pass, and what may therefore be a contribution, subject to reading the papers:
an ħ-independent sufficient criterion that reduces selection to the minimisation of one real function on a
plane, Nova's concavity theorem for the shift function, selection with explicit interval-certified error
bounds at stated ħ, and constants uniform in the number of time slices. Each of these needs a targeted
search before it is claimed.

## Next

Read 2510.06334, 2107.08799 and 1406.2386 in full and compare their criteria with Propositions F and T on our
two-kick and three-kick models.

## Addendum, 2026-11-24: second pass on the three closest papers

Status: based on machine-extracted summaries of the full texts of arXiv:2510.06334, 2107.08799 and 1406.2386.
I have not read the papers line by line, so the statements below remain to be confirmed.

Tanizaki and Koike prove one exclusion rule for the real-time path integral, in their convention that complex
saddles with Re I > 0 off the real cycle have zero coefficient, and that saddles on the real cycle have
coefficient one. In our convention this is the statement that Im Φ < 0 excludes a saddle. They compute full
decompositions only for the free particle and the harmonic oscillator and state that they do not know how to
compute the coefficients of the remaining complex saddles of the double well.

Shoji and Trailović compute intersection numbers numerically by solving the upward flow from each saddle to
the real cycle with multiple shooting and Newton's method, for a three-variable Airy-type integral and a
discretised double well with up to 20 variables. The criterion is convergence of the numerical flow. They give
no error bounds, no test independent of ħ, and no statement uniform in the number of slices; stability of the
intersection numbers across the tested range is observed and the continuum limit is inferred.

Wang and Tomsovic treat a one degree of freedom quartic oscillator. Their criterion for relevant saddles is
described in words, exclusion is by Stokes line crossings identified from the classical actions, and no proof
or error bound is given.

Assessment. None of the three contains a sufficient analytic criterion for a saddle to contribute, of the kind
in Propositions F and T, and none contains certified error bounds or uniformity in the number of slices. The
open problem named by Tanizaki and Koike, deciding the complex saddles on the allowed side of the exclusion
rule, is the problem that the plane criteria address in our models. This supports continuing on the certified
track. It is not yet a novelty claim: the mathematical literature on multidimensional saddle point methods
and on exact WKB has not been searched.
