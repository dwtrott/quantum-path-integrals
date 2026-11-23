# Research correspondence 007: Sparks → Nova and David
**Date:** November 23, 2026
**Project:** Certified Trajectory Selection in Quantum Path Integrals
**Subject:** Proposal for the multi-saddle direction, with a scoping experiment
**Status:** Proposal and numerical evidence. No theorem, no proof, no Lean. Prior art not yet checked.

## 0. Summary

David prefers the more novel direction over an extension of Hargé's single-saddle result. This note defines the first multi-saddle problem concretely, reports a scoping experiment, and proposes a route from one variable to $N$ slices.

The first configuration with several real saddles in our model is reached with an attractive Gaussian well and is organized by a cusp caustic in the coupling and endpoint parameters. Outside its fold line there is one real saddle, one complex saddle that must be kept, and infinitely many complex saddles that must be bounded. That is the trajectory selection problem of the original handoff in its smallest exact form.

## 1. The one-variable problem

Take $N=2$, so one interior variable, and write the endpoint-stripped integral as
$$I(\hbar)=\lim_{\eta\to0^+}\int_{\mathbb R}e^{i\varphi(x)/\hbar-\eta x^2}dx,\qquad\varphi(x)=\tfrac12\lambda x^2-cx-g\,e^{-x^2/2\sigma^2},$$
with $\lambda=(m/d)(2-\omega^2d^2)$, $c=(m/d)(a+b)$ and $g=\epsilon dV_0$. A negative $g$ is an attractive well.

**Exact reference.** The Fourier series of correspondence 002 gives $I$ exactly as an entire series in $g$,
$$I=\sum_{p\ge0}\frac{(-ig/\hbar)^p}{p!}\sqrt{\frac{2\pi}{p/\sigma^2-i\lambda/\hbar}}\exp\!\Big(-\frac{(c/\hbar)^2}{2(p/\sigma^2-i\lambda/\hbar)}\Big),$$
which can be evaluated to any precision and later enclosed in interval arithmetic. The multi-saddle regime lies far outside the tree domain, so the series is a reference and not an analytic tool there.

**Saddle structure.** At $c=0$ the origin is a saddle with $\varphi''(0)=\lambda+g/\sigma^2$. For $\lambda>0$ it changes sign at $g=-\lambda\sigma^2$, where two further real saddles appear. With $c$ as the second parameter this is a cusp. Inside the cusp region there are three real saddles. Across the fold line two of them merge and continue as a complex conjugate pair. In addition, because $\varphi$ is entire with an essential singularity at infinity, there are infinitely many complex saddles, near $x^2\approx-2\sigma^2(\log(-\lambda\sigma^2/g)+2\pi in)$.

## 2. Scoping experiment

Parameters $\lambda=\sigma=1$, $g=-2$. The fold point is $x_f=0.62594$ and the fold line is at $|c|=0.40322$. Each saddle contributes its leading term, with the descent direction oriented left to right. Reference values from the series at 120 digits. Script: `experiments/2026-11-23-sparks/two_saddle_toy.py`.

**Outside the fold**, $c=|c_{\rm fold}|+0.15$: one real saddle at $x=1.515$; the nearest complex pair at $-0.634\pm0.327i$ with $\operatorname{Im}\varphi=\pm0.0330$; the next pairs have $|\operatorname{Im}\varphi|\approx5.0$ and $7.6$.

| $\hbar$ | $|I|$ | Error, real saddle only | Error, real saddle plus the decaying complex saddle | $e^{-\operatorname{Im}\varphi/\hbar}$ |
| --- | --- | --- | --- | --- |
| 0.1 | 0.980 | 0.514 | 0.0701 | 0.719 |
| 0.05 | 0.594 | 0.274 | 0.0229 | 0.517 |
| 0.025 | 0.330 | 0.103 | 0.00545 | 0.267 |
| 0.0125 | 0.202 | 0.0197 | 0.00084 | 0.0713 |

Keeping only the real saddle leaves an error of 10 to 50 percent of $|I|$ at these values of $\hbar$. Adding the single complex saddle with $\operatorname{Im}\varphi>0$ reduces it by a factor of 7 to 23. The conjugate saddle, with $\operatorname{Im}\varphi<0$, would contribute an exponentially large term and is correctly excluded. The same pattern holds at $c=|c_{\rm fold}|+0.05$, where the improvement factor is 3 to 12 because the saddles are closer to coalescence.

**Inside the cusp region**, $c=0.2$: three real saddles at $-1.013$, $-0.209$, $1.311$. The best single saddle leaves errors of 0.11 to 0.65. The three-saddle sum leaves 0.28, 0.12, 0.0094, 0.026 at $\hbar=0.2,0.1,0.05,0.025$. The decrease is not monotone, which I attribute to interference among the first corrections of the three saddles and to the proximity of two of them. This needs the next-order terms before any rate is stated.

## 3. Target statements

1. **Selection outside the fold.** For $(c,g)$ at distance at least $\delta$ from the fold line, $|I-I_{\rm real}-I_{z_1}|\le B(\hbar)$ with explicit $B$, where $z_1$ is the adjacent complex saddle with $\operatorname{Im}\varphi>0$ and $B$ consists of an explicit bound of order $\hbar^{3/2}$ on the retained expansions plus a certified bound on the sum over all other complex saddles. The statement includes a proof that $z_1$ has coefficient one and its conjugate coefficient zero.
2. **Three real saddles.** Inside the cusp region and away from its boundary, the three-saddle sum with an explicit remainder.
3. **Uniform statement through the fold.** An Airy-type uniform approximation with explicit error. OscInt already contains fold and cusp decay bounds with explicit constants, which may serve the Lean side.

Each statement concerns the finite-dimensional integral with the regulator removed.

## 4. From one variable to N slices

For general $N$ the onset of multiple real saddles occurs when the smallest eigenvalue in modulus of the perturbed Hessian $A_\epsilon$ passes through zero. Near that point one direction is soft and $N-2$ directions are stiff. The proposal is a certified reduction:

1. Fix the soft coordinate $u$ and integrate over the stiff directions. Conditional on $u$ the stiff Hessian is uniformly nondegenerate, so the single-saddle machinery of correspondence 003 to 005 applies with $u$ as a parameter and yields an effective phase $\varphi_{\rm eff}(u)$ and amplitude with an explicit relative error, uniform in $N$.
2. Apply the one-variable statements of Section 3 to the effective integral.

This is the real-time counterpart of the reduction to quasi-zero modes in Dersy and Schwartz, and it gives the single-saddle work a role as a component instead of an end result.

## 5. Prior art to check before any claim

In one variable, error bounds for steepest descent are an established subject. The audit must cover Berry and Howls on hyperasymptotics for integrals with saddles, Boyd's error bounds for the method of steepest descents, the explicit bounds of Nemes and of Bennett, Howls, Nemes and Olde Daalhuis, and uniform approximations at coalescing saddles (Chester, Friedman, Ursell). I cite these from memory and have verified none of them in this session. On the path integral side: Feldbrugge's resurgence paper, and Picard–Lefschetz treatments of real-time path integrals. A possible gap is the combination of an entire phase with infinitely many saddles and a certified tail, the reduction of Section 4 with constants uniform in $N$, and formal verification. A one-variable result alone is unlikely to be new.

## 6. Proposed first tasks

- Sparks: next-order corrections for the three saddles and a first explicit bound for target 1 in one variable, with the saddle coefficients proved by an explicit contour.
- Nova: independent derivation of the saddle structure and coefficients for the same phase, and the prior-art audit of Section 5. Nova's branch `nova/three-real-saddles` suggests work in this direction has begun; we should compare only after both derivations exist.
- David: choose between the Gaussian well and the Rosen–Morse profile as the multi-saddle model. The Gaussian has an exact series and infinitely many saddles; Rosen–Morse has poles, an exact propagator and is the benchmark in Feldbrugge.

— Sparks
