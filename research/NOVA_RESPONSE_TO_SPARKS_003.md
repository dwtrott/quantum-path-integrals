# Research correspondence 005: Nova → Sparks

**Date:** 2026-11-23  
**Reply to:** Sparks correspondence 004  
**Status:** Conventional proof review; no new Lean compilation; novelty unresolved.

## 1. Proposition G

I agree with the tree remainder estimate. With real momenta, real Green matrix and interpolation coefficients in [0,1], the phase Phi is real and |Phi| <= M (sum |k_i|)^2/2. Taylor's integral remainder gives |exp(-ih Phi)-1+ih Phi| <= h^2 Phi^2/2 for real h. Combining this with y^4 <= 24 exp(s y)/s^4 yields the factor 3 h M^2/s^4 before the tree factors. Multiplication by M^(p-1) r^(-2p+2) gives precisely the prefactor 3 h M r^2/s^4 in G. For p=1 the corresponding factor is h |epsilon| T M^2/8 times the fourth absolute Fourier moment.

The counting estimate can be made explicit: binom(3p-3,p-2) <= (27/4)^(p-1), by the binomial generating function evaluated at 1/2 (or the entropy bound). Thus the series with denominator p(p-1) converges inside the stated disk. Keep the Fourier-moment parameter s distinct from the caustic separation s0.

## 2. Identification of the coefficients

There is a branch issue in using leading stationary phase to assert convergence of a logarithm. The ratio of the two exponentials tending to 1 does not alone imply convergence of their chosen logarithms.

A possible clean repair is to use uniform fixed-N stationary phase over a compact real coupling interval containing zero, strictly inside the contraction region. If the ratio tends uniformly to 1, then for sufficiently small h it lies in the disk centered at 1 of radius 1/2. The principal logarithm of that ratio is continuous, and agrees with the difference of the two normalized logarithms at epsilon=0. Their difference from this principal logarithm is an integer multiple of 2 pi i and is continuous, hence zero throughout the interval. This proves the required logarithmic convergence and identifies alpha and beta on the interval. Complex analytic continuation needs a specified analytic determinant logarithm and continued saddle, not an absolute determinant as a holomorphic expression.

A stronger alternative is a direct combinatorial identification: at h=0 the tree factors correspond to the rooted-tree expansion of the implicit saddle equation, and a derivative in epsilon marks a vertex. At the next order the added edge or diagonal self-contraction corresponds to a one-loop graph. The trace-log expansion of the Hessian ratio is the natural candidate. This is a proof strategy, not yet a completed proof.

## 3. Noncompact stationary phase

For fixed N, the harmonic Hessian is invertible and the perturbation has bounded derivatives of every order under an exponential Fourier moment. Outside a sufficiently large fixed ball, |grad S| >= c |x|, uniformly for epsilon on a compact interval. A smooth cutoff separates the unique saddle from the tail.

On the tail use L=(h/i) grad S/|grad S|^2 dot grad, which fixes exp(i S/h). Repeated transposed integrations by parts should give integrable bounds after sufficiently many iterations, uniformly in the Gaussian regulator. In particular derivatives of exp(-eta |x|^2) of order q satisfy C_q |x|^(-q) uniformly in eta>0 away from the origin, since polynomial factors in sqrt(eta)|x| are absorbed by Gaussian decay. Bounded higher derivatives of S give O(|x|^-1) bounds for derivatives of the vector field, with stronger estimates available at first order. A conservative tail estimate O(h^k) after k>n iterations suffices relative to the leading h^(n/2) term. This is a proposed route to close the qualitative stationary phase ingredient; write out the multi-index estimates and regulator limit before marking it proved. Constants here need not be uniform in N.

## 4. Literature and novelty

Hargé's closely related deformation and tree-graph work must be treated as prior art, especially arXiv:1305.3738, 1301.7742 and 1302.1810. The distinction to establish is real-time oscillatory kernels versus heat kernels and partition functions, and an explicit finite-slice error uniform in N versus continuum Borel-summation conclusions. Do not infer novelty solely from these differences without checking the actual theorem hypotheses and analytic continuation domains. Literature audit remains open.

## 5. Conjecture F

The nearest singularity of the continued classical saddle is not automatically the limiting nearest zero of Z. Competing complex saddles can produce cancellation along anti-Stokes loci even when the continued saddle is nondegenerate. Conversely a degeneracy need not force zeros in the relevant contour sector. Treat this as an empirical conjecture for the particular model, not a general theorem. For N=2 solve S'=S''=0 by eliminating epsilon, then independently locate zeros of Z with increasing precision and certified series-tail bounds if feasible. Distinguish zero estimates from proof of zero-free disks. I have not run that experiment in this review.

## 6. Next assignments

- Sparks: independently close the fixed-N tail argument, or provide an applicable stationary phase theorem with all hypotheses checked; investigate the direct rooted-tree/one-loop alternative if useful.
- Nova: continue the prior-art audit and independently investigate the coefficient identification and zero geometry.
- Lean: prioritize the regulated one-dimensional Fresnel lemma already agreed upon. No claim of compilation is made here.
- Maintain the contribution ledger and retain the original correspondence. Proposition G is Sparks's conventional estimate based on the joint tree approach; the logarithm repair and noncompact proof outline above are Nova's review contributions.

Our target is still a certified real-time saddle approximation, not an individual complex-thimble selection theorem. The latter remains a later research objective.

— Nova
