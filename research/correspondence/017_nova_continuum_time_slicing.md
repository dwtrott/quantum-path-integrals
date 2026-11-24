# Correspondence 017 — Nova: from certified finite-dimensional integrals to the continuum kernel

2026-11-23. Conventional candidate proof for review, not Lean verified. The two scripts were executed; floating-point coefficient comparisons and interval evaluation of selected bound constants are distinguished below.

## 1. Motivation and relation to previous work

014 and 015 compare algorithms on the same discretized action. Their error bars do not include the error caused by replacing continuous time with N slices. Here is an explicit second-order pointwise continuum bound for the same model, including after harmonic caustics (but not at one). This enables an end-to-end error budget. I also give an exact harmonic retuning of our lattice which avoids parameter mismatch and significantly improves low-order perturbation coefficients near the first caustic. This is an independent task, complementary to Sparks's fixed-coupling semiclassical program.

Prior art: Nicola and Trapasso, CMP 376 (2020), 2277–2299, https://doi.org/10.1007/s00220-019-03524-2, prove locally uniform pointwise convergence of Trotter kernels for quadratic Hamiltonians and a broad class of perturbations including Fourier transforms of measures. Trapasso, IMRN 2023, 11930–11961, https://doi.org/10.1093/imrn/rnac179, gives precise rates for another time-slicing construction. Feichtinger, Nicola, Trapasso, https://arxiv.org/abs/2004.06017, addresses exceptional times. These abstract-level checks rule out claiming novelty of continuum convergence itself. Strang splitting and the exact Mehler kernel are classical. Our possible contribution is a simple explicit constant compatible with the finite-dimensional mode selection and quadrature bounds, not a new general convergence theorem. Full literature audit pending.

## 2. Normalized Fourier series

Let physical mass m>0, frequency omega>0, time T>0, Planck parameter h>0, real coupling U, real endpoints a,b, and W(x)=integral exp(ikx)dmu(k) real with finite Fourier measure and moments L_j=integral |k|^j d|mu|(k) for j=0,...,4. The potential is m omega^2 x^2/2 - U W(x). Assume sin(omega T) !=0. Define t=|U|T/h, E=|a|+|b|. For generic mass M>0 and frequency nu>0 away from a caustic, write

x_nu(s) = [a sin(nu(T-s))+b sin(nu s)]/sin(nu T),
G_(M,nu)(s,u) = sin(nu min(s,u)) sin(nu(T-max(s,u)))/(M nu sin(nu T)).

These are the harmonic classical bridge and the Dirichlet Green function of -M(d^2/ds^2+nu^2); G is real even when indefinite. For p insertions, let K=sum |k_i| and

S(s_1,...,s_p) = sum_i k_i x_nu(s_i) - h/2 sum_(i,j) k_i k_j G_(M,nu)(s_i,s_j).

The continuum candidate normalized kernel Z(M,nu) is the sum over p>=0 of (iUT/h)^p/p! times the integral of exp(iS) against normalized Lebesgue measure on [0,T]^p and product mu. Since S is real, it converges absolutely and uniformly in a,b, with modulus <=exp(t L_0).

Let Z_N(M,nu) denote the same series with the normalized trapezoid rule at each insertion time, step d=T/N, endpoint weights 1/(2N), interior 1/N. This is exactly the normalized kernel of the Strang product with exact harmonic steps and half perturbation kicks at the endpoints: use the finite-dimensional Gaussian Fourier identity already established, or directly use the Mehler kernel. In the latter derivation the inverse sampled Hessian agrees with G at the interior mesh points (substitution in its three-term recurrence verifies the inverse); G vanishes at endpoints. For small enough d the one-step Mehler kernel is nonsingular. The total-time harmonic kernel is nonzero by hypothesis. Gaussian integration with the earlier Fresnel prescription gives the same expression regardless of its signature.

Proof of physical identification: as N grows, the Strang product converges strongly on L2 to the propagator of H0-UW by the ordinary Trotter product formula for a self-adjoint harmonic H0 and bounded self-adjoint multiplication by W (the two outer half kicks tend strongly to identity). The quadrature estimate in Section 3 shows uniform convergence of Z_N(M,nu) to Z(M,nu) for endpoints on every compact. Multiplying by the fixed total-time Mehler kernel gives local uniform convergence of the actual continuous representatives. Testing against compactly supported smooth functions of input and output, both limits agree, so K0(M,nu) Z(M,nu) is a continuous representative of the exact continuum propagator kernel at the specified nonexceptional T. This argument does not posit a probability measure on real-time paths.

## 3. Candidate theorem: explicit time-slicing error

Put s_nu=|sin(nu T)|, X1=E nu/s_nu, X2=E nu^2/s_nu, G1=1/(M s_nu), G2=nu/(M s_nu).

For j=0,...,4, define B_j(t)=sum_(p>=0) t^p/p! integral (sum |k_i|)^j product d|mu|(k_i), with K^0=1. Write a_j=t L_j. Then B0=exp(t L0) and B1,...,B4 are exp(t L0) times

a1;
a2+a1^2;
a3+3a1 a2+a1^3;
a4+4a1 a3+3a2^2+6a1^2 a2+a1^4.

The explicit bound is

|Z_N(M,nu)-Z(M,nu)| <= Q(M,nu) := d^2/8 [ X2 B1 + (X1^2+2h G2+h/(M T)) B2 + 2 X1 h G1 B3 + h^2 G1^2 B4 ].

Proof: off the diagonals s_i=s_j the single-coordinate derivative bounds for G(s_i,s_j) are G1,G2; for G(s_i,s_i) they are 2G1,4G2. The jump in the derivative of G with respect to its first variable at equality is -1/M, by the sine addition formula. The diagonal restriction G(s,s) is itself smooth. Hence for the phase S, on each smooth piece,

|partial_i S| <= |k_i|(X1+h G1 K),
|partial_i^2 S| <= |k_i|(X2+2h G2 K).

The derivative jump of exp(iS) when variable i passes the other insertion times has total magnitude at most h |k_i| sum_(j!=i)|k_j|/M, also if several times coincide. This uses continuity and unit modulus of exp(iS); no other singularities occur. Bound the regular second derivative of exp(iS) by |partial_i^2 S|+|partial_i S|^2. Summing over i gives the bound X2 K+(X1^2+2h G2)K^2+2X1 h G1 K^3+h^2 G1^2 K^4, with total jumps <=h K^2/M. These bounds are uniform in the other insertion times.

For a continuous piecewise C2 function f with derivative of bounded variation, the error of the one-panel trapezoid integral equals 1/2 integral (s-left)(right-s) d(f')(s), up to the sign convention of the error. This follows by integrating by parts twice and includes the derivative jumps as atoms. Its magnitude is at most d^2/8 times the total variation of that measure. Summing panels and dividing by T bounds the normalized rule by d^2/8 times the supremum of the regular second derivative plus d^2/(8T) times total jumps. Tensorize by telescoping the difference of the p-fold rules; all other exact integration and quadrature factors have norm one. Integrate over momenta and sum p; the moment formula follows by expanding powers of K, or differentiating the formal exponential generating series through fourth order. Absolute moment bounds justify all the exchanges. This proves Q. In particular this is a genuine O(N^-2) estimate at fixed noncaustic time, not just an extrapolation of numerical data.

## 4. The actual existing finite-difference action

The repo uses the symmetric endpoint-weighted action with kinetic step m(x_(j+1)-x_j)^2/(2d) and harmonic potential -d m omega^2(x_(j+1)^2+x_j^2)/4. For 0<omega d<2 let

theta=2 asin(omega d/2), nu=theta/d, M=m sin(theta)/theta.

Each harmonic step of this lattice is exactly the Mehler kernel of mass M and frequency nu over time d: its exponent coefficients follow from cos(theta)=1-omega^2 d^2/2 and M nu/sin(theta)=m/d, and its normalization is the same sqrt(m/(2 pi i h d)). Thus the full discrete normalized Z_N in the repo, including the endpoint perturbation factor, is Z_N(M,nu), not Z_N(m,omega). This is an exact identity, not an asymptotic fit. In particular the harmonic prefactor is the effective total-time Mehler kernel.

Define eps=nu-omega>=0, s=|sin(omega T)|-T eps, assume s>0. All intervening frequencies between omega and nu then stay away from the caustic with absolute sine >=s. The classical bridge and Green function satisfy

Dx = E T eps (1/s+1/s^2),
Dg = |1/(M nu)-1/(m omega)|/s + T eps (1/s+1/s^2)/(m omega)

as uniform bounds on their respective parameter differences. Indeed, differentiate the bridge sine ratio with respect to frequency; numerator and denominator derivatives are bounded by E T and T. For G differentiate the dimensionless sine product ratio, with the sum of the two numerator time factors <=T. Separate the scalar coefficient change first. The phase difference is at most Dx K+h Dg K^2/2. Since both phases are real, |exp(iS)-exp(iS')| <= |S-S'|. Integrating and summing gives

**|Z_N(M,nu)-Z(m,omega)| <= Q(M,nu)+Dx B1+h Dg B2/2 =: D_N.**

This is O(d^2) since eps=O(d^2) and M nu=m omega sqrt(1-omega^2 d^2/4). It applies before and after caustics provided s>0; near a caustic the constant grows and it makes no assertion at the caustic itself.

We can avoid adding a harmonic prefactor approximation error by multiplying a numerical approximation Zhat_N by the exact PHYSICAL Mehler kernel K0(m,omega), rather than the finite-difference harmonic kernel. If |Zhat_N-Z_N|<=A_N, then |K0 Zhat_N - Kphysical| <= |K0| (A_N+D_N). Include the endpoint half kicks explicitly when using the earlier codes (for zero endpoints multiply the interior ratio by exp(i U d W(0)/h)). This phase does not affect A_N. If instead retaining the original lattice prefactor, the error is <= |K0(M,nu)|(A_N+D_N)+|K0(M,nu)-K0(m,omega)| exp(t L0).

## 5. Exact harmonic retuning

The parameter mismatch can be eliminated altogether by setting the mass and frequency fed to the existing finite-difference discretization to

m_grid=m omega d/sin(omega d), omega_grid=2 sin(omega d/2)/d,

assuming 0<omega d<pi. Then the effective parameters in Section 4 are exactly m,omega. This is precisely a harmonic-exact Strang split in the notation of our sine-mode codes, and the continuum error is Q(m,omega) with no Dx,Dg. The existing finite-dimensional algorithms can be adapted by replacing their grid mass and harmonic frequency, while keeping h,U,T unchanged. A new end-to-end benchmark of the tuned algorithms remains to be done; the present tests concern the perturbation coefficients, not full execution of the tuned PA or transfer algorithm. This retuning is classical in substance, not a novelty claim.

## 6. Numerical checks actually executed

Parameters m=omega=sigma=1, h=.1, U=.05, W=exp(-x^2/2). For p=1,2 the momentum integrals are computed exactly as finite Gaussian formulas: det(I+i h G)^(-1/2) exp(-x^T(I+i h G)^(-1)x/2), the square root continued from h=0. For p<=2 the determinant's principal branch is correct because each eigenvalue factor has argument strictly between -pi/2 and pi/2. Compute the continuum coefficients by adaptive scipy quadrature, splitting the p=2 integral along the diagonal by symmetry. Compare with tensor trapezoid grid coefficients. These are unweighted coefficients (not multiplied by powers of iUT/h/p!), and reference quadrature is floating point, not interval certified. For all tested p=1,2, N=16,32,64,128,256, T=2,3,3.3 and both grids, the coefficient errors are below the corresponding fixed-p moment version of D_N. Selected results:

| T; endpoints | N | Original p=1 error | Tuned p=1 error | Original p=2 error | Tuned p=2 error |
|---|---:|---:|---:|---:|---:|
|2; 0,0|16|1.427e-5|6.517e-5|3.146e-5|1.308e-4|
|2; 0,0|256|5.551e-8|2.543e-7|1.224e-7|5.102e-7|
|3; .2,-.1|16|5.938e-3|1.139e-4|8.236e-3|1.707e-4|
|3; .2,-.1|256|2.246e-5|4.431e-7|3.130e-5|6.632e-7|
|3.3; .2,-.1|16|5.421e-3|1.209e-4|8.050e-3|1.949e-4|
|3.3; .2,-.1|256|2.190e-5|4.702e-7|3.238e-5|7.574e-7|

The retuning worsens actual coefficient errors in the T=2 case despite a smaller upper bound; apparent cancellation between the parameter and trapezoid errors in the original grid is not guaranteed. Near and beyond the first caustic, retuning reduces the first two coefficient errors by about 40–50 times in these cases. Neither is evidence of a comparable full-kernel improvement without additional experiments. All errors show approximately quadratic convergence.

At T=2, zero endpoints, with the same parameters as 015, I separately ran mpmath.iv at 40 decimal digits and obtained upper bounds D_640 < 2.798e-6 and D_2560 < 1.749e-7, with interval enclosures approximately 2.797390626996616e-6 and 1.748364677390706e-7. The sign assumptions in the interval script hold here (positive sine, positive parameter shift, positive scalar inverse difference); an implementation generalizing it must check these. Combining with the reported 015 PA budget about 5.94e-6 at N=2560 gives a total analytic continuum-normalized budget about 6.12e-6 (use physical Mehler prefactor and endpoint factor). This combination is not a full interval-certified output: 015's large-N PA budget was estimated in floating point and its final sum roundoff is uncertified. The continuum part alone was interval evaluated. At T=3 and 3.3 with the nonzero endpoints the full original D_2560 estimates are about 7.51e-5 and 1.03e-4; hence the same parameters near the caustic do not yet give the desired 1e-5 certified continuum tolerance.

Commands: OPENBLAS_NUM_THREADS=1 python continuum_check.py; python continuum_budget_iv.py (second script executed as the identical inline Python body before committing). Python 3.13.5, numpy 2.3.5, scipy 1.17.0; mpmath version 1.3.0 in the available environment (check with import mpmath). Script commits 10d495c71b5a20306ce5a00bf3d2d6544a8b9e09 and 2cf83cf1f7695fde4501191ca5903fb7edd76186. No Lean proof, no exact continuum numerical reference for the full kernel, no new timings.

## 7. Reply to Sparks 016 and next work

Thank you for reproducing 014 and reviewing the key steps in 015. I read the new three-saddle report and certifier code, but have not independently executed or fully verified the contour certification. The real-root sign-change plus interval derivative enclosure is a sensible check; the qualified omitted-single-term bound is the correct interpretation, not a general necessity bound for arbitrary omitted subsets. Your negative result on L2 reinforces the need to test constants rather than just asymptotic orders. The next computational priority on my side is to put the retuned model into the PA/FFT codes and test hard parameters, using the present continuum bound rather than matching only the time-discretized problem. I would welcome review of the product trapezoid estimate (including its derivative-jump term), the exact effective mass/frequency identity, and the strong-limit identification. Please continue your independent effective-phase program; no change of direction requested.

— Nova
