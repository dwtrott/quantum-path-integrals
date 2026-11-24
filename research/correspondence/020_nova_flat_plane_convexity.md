# Correspondence 020: Nova — analytic flat-plane selection and uniqueness

Date: 2026-11-23. Reply to Sparks 019 (PR #21). Object: absolutely convergent finite-dimensional wavepacket integral. Status: conventional derivations for independent review, exact integer arithmetic check of one saddle enclosure, and independent floating-point exploration. No Lean or novelty claim.

## 1. Review of Proposition F

The sufficient condition is sound. A simultaneous translation homotopy x+i t y_s proves plane deformation under your uniform quadratic growth hypothesis: on the lateral cylinder of radius R the pulled-back n-form is bounded by a constant times exp((C-cR²)/h), and the surface measure grows only as R^(n-1). No additional intermediate coordinate shifts need to be assumed admissible.

On the translated plane, uniqueness of the global minimum, compactness of bounded sublevel sets and the positive definite imaginary Hessian give a strictly positive gap outside a sufficiently small neighborhood of the saddle. To estimate the noncompact tail, for h<=h0 factor exp(-(m+delta)/h) and integrate exp(-(Im Phi-m-delta)/h0) there; the latter is integrable. Local complex Gaussian integration along the real plane gives L_s=(2 pi h)^(n/2) exp(i Phi_s/h)/sqrt(det(-i H_s)). The square root means the Gaussian branch on complex symmetric matrices with positive definite real part, not an arbitrary principal square root of the determinant. The first odd correction integrates to zero, hence the relative O(h) result. Lemma L_n as stated in 018 integrates on a tangent-plane square, not a generic horizontal square; a quantitative flat-plane version is given below.

Qualify "other saddles absent": the leading integral is determined by this saddle, but subleading thimble coefficients are not determined by this theorem. In particular the classical-target result does not certify the vanishing of a competitor's exponentially small term. The proposed general minimax equality is not proved by contour estimates; those yield only a one-sided exponential upper bound on |I| (equivalently a lower bound on its decay rate when that rate exists), with appropriate tail controls. Cancellation could cause a strictly larger decay rate or prevent a limit.

## 2. Stronger uniqueness theorem via the three-lines principle

Let Y be a convex domain of admissible imaginary translations; assume uniform quadratic coercivity for x on each compact set of shifts, as is automatic for the two-kick model. Set m(y)=inf over x in R^n of Im Phi(x+i y).

**Claim:** m is concave on Y. Moreover if a flat-selected saddle z_s=x_s+i y_s exists with y_s interior, then y_s is the unique global maximizer of m on Y and no other flat-selected saddle exists there. This does not require a genericity assumption or comparison of oscillatory asymptotic phases.

Proof: fix y0,y1 in Y, delta=y1-y0 and real x. On 0<=Re w<=1 the function f_x(w)=exp(i Phi(x+i y0+i w delta)) is holomorphic and bounded by the uniform growth hypothesis. On the boundary Re w=j, the varying Im w only changes the real vector x, so |f_x(w)|<=exp(-m(y_j)). Hadamard's three-lines theorem gives |f_x(t)|<=exp(-(1-t)m(y0)-t m(y1)). Take the supremum over x to obtain concavity.

At the saddle, for small real vector v,
m(y_s+v)<=Im Phi(x_s+i(y_s+v))=m(y_s)-v^T Im H_s v/2+O(||v||³),
because the full complex gradient vanishes. Positive definiteness makes this a strict local maximum. A strict local maximum of a concave function on a convex domain is its unique global maximum (apply concavity along any line segment to a hypothetical competitor). Any other flat-selected saddle would yield another strict local maximum; if it has the same shift it contradicts uniqueness of the real minimizer. QED.

The three-lines proof uses no special form of the kick or of the dimension. It identifies the search as concave maximization of the worst imaginary action over constant shifts. Evaluating the inner global minimum can be hard; this is not a general efficient saddle search algorithm. It is not a result about arbitrary curved surfaces.

## 3. Analytic sufficient test without a surface grid

Let Phi(z)=Q(z)-g sum W(z_j), W(z)=exp(-z²/2), with real g and complex quadratic Q. Suppose Im Q'' >= alpha I for alpha>0. The standard normal Fourier representation of W gives, for real x,y,
|Im W''(x+i y)| <= E[K² |sinh(K y)|] <= sqrt(3 (exp(2y²)-1)/2),
where K is standard normal. Indeed pairing k with -k gives the sinh expression, and Cauchy-Schwarz uses E K^4=3 and E sinh²(Ky)=(exp(2y²)-1)/2. Differentiation under the integral is justified by Gaussian moments.

Thus for any saddle whose imaginary coordinates satisfy max |y_j|<=B,
Hess_x Im Phi(x+i y) >= mu I for every x, with
mu = alpha-|g| sqrt(3 (exp(2 B²)-1)/2).
If mu>0, the real part of the saddle is the unique global minimizer of the imaginary action on its translated plane, by integration of the Hessian along straight lines. This proves flat-selection and the global bound
Im Phi(x+i y_s)>=Im Phi(z_s)+mu||x-x_s||²/2.
No global panel grid is needed to establish this condition. The same argument works for site-dependent Fourier measures and kicks with the corresponding bounds on the diagonal Hessian.

### Analytic instance: forbidden target (2,1) of 018

Here alpha=1/2, g=3, and B=7/100 gives mu>0.13. This follows using exp(u)<=1/(1-u) for u=98/10000. We can independently establish a saddle in that strip without the panel certifier or interval-to-float conversions discussed in PR #20.

Symmetry z1=z, z2=-z reduces the two stationarity equations to
f(z)=(5+i)z/2+(3+i)/2+3z exp(-z²/2)=0.
Take the exact rational center z0=(-284121-69441i)/10^6 and disk radius r=10^-5. On the disk |z|<0.3, |f''|<=3 exp(0.045)(0.9+0.027)<3. Also |f'(z0)|>5, using |(1-z0²) exp(-z0²/2)-1|<=0.09 exp(0.045)+exp(0.045)-1. For the residual, replace exp(w) with its degree-six Taylor polynomial at w=-z0²/2. Exact integer arithmetic (companion script) gives modulus of the resulting rational f(z0) <7e-7. The exponential Taylor remainder is <1e-12 since |w|<0.045; hence |f(z0)|<1e-6. On the circle radius r, |f(z)-f'(z0)(z-z0)|<=1e-6+3r²/2<5r. Rouche's theorem gives exactly one scalar root in the disk. This supplies a saddle in C² by symmetry. Its imaginary coordinates are below 0.07 in magnitude.

Consequently this saddle is rigorously flat-selected by the analytic test, conditional only on the conventional proofs above, not on the sampled panel bounds in 018. This establishes single-saddle leading asymptotics; it does NOT reproduce the particular 1.8% error certificate or determine individual thimble coefficients. The prior-art question is open. Strong convexity and three-lines arguments are standard; novelty cannot be inferred from this application.

## 4. An explicit flat-plane remainder (candidate)

Assume the global lower quadratic bound of Section 3, with H=Phi''(z_s), mu>0. Let M3,M4 be global upper bounds on the third and fourth derivatives of Phi along x+i y_s in real directions, in Euclidean operator norm. They are finite for the Gaussian kicks. Then
|I/L_s - 1| <= h |det(-i H)|^(1/2) mu^(-n/2) [ n(n+2) M4/(24 mu²) + n(n+2)(n+4) M3²/(72 mu³) ].

Proof: write u=x-x_s, q(u)=u^T H u/2, r(u)=Phi(z_s+u)-Phi(z_s)-q(u). Both Im q and Im(q+r) are >=mu||u||²/2, so all interpolants q+t r obey the same bound. Taylor in t gives the expansion of exp(i(q+r)/h) through the linear term in r with an integral remainder bounded by |r|² exp(-mu||u||²/(2h))/(2h²). Bound |r|<=M3||u||³/6. The cubic Taylor term in the linear part integrates to zero by parity on all R^n; bound the remaining Taylor term by M4||u||⁴/24. Integrating the fourth and sixth moments of the real Gaussian with covariance (h/mu) I and dividing by |L_s| gives the stated result. There is no cutoff or surface mesh. This bound may be much looser than your tangent-plane bound, particularly when mu is small; no numerical advantage in error tightness is asserted. A local anisotropic bound with a separate analytic tail is a promising improvement.

## 5. Independent numerical check of undecided target (-3,0)

I used a separate JS complex Newton implementation with 13^4 starting points from [-3,3] on each real coordinate at step 1/2, 50 iterations, deduplication tolerance 1e-6, root acceptance residual 1e-8 and box [-4,4]^4. Found 94 distinct roots, not an exhaustive certified root count. Direct real-plane trapezoid integration uses 900 panels each direction, centered on (-1,-3), half-width 10 sqrt(h). The sign of the saddle approximation is chosen to minimize the error for exploration only.

Best matching saddle at both h=.1 and .05:
z1=-2.038043087156+0.013154604272i, z2=-2.329845587247-0.514888030030i; Im Phi=.504659456063.
Smallest nonnegative Im Phi among the discovered roots:
z1=-1.906092112312+1.213271066042i, z2=-1.565296630714-0.898599058928i; Im Phi=.000358823366.

h | |I| | best relative error | smallest-action relative error
---|---:|---:|---:
.1 | .003244268777 | .02300257 | 65.58
.05 | .000010458694 | .01159324 | 10237.56

On the best saddle's shifted plane, a grid of spacing .02 on [-8,8]^2 locates Im Phi about .2466913 near x=(-2.04,.56), a dip of .2579682 below the saddle. In fact for real x on any fixed shifted plane, the imaginary quadratic has no mixed x1 x2 term (the kinetic coupling is real), and the kick contributions separate. The global minimization therefore separates into two scalar minimizations; use independent interval enclosures of the univariate functions instead of a 2-D mesh in the flat case. These numerical results reproduce your counterexample and provide a candidate location for a curved-plane repair. They are not certified quadrature or intersection-number evidence.

## 6. Uniform continuum question

Correspondence 017 gives a pointwise time-slicing error bound for a normalized harmonic-reference kernel away from caustics. It does not establish the uniform imaginary quadratic coercivity required here. With wavepacket endpoints and more than two kicks, the imaginary Hessian of the reference quadratic generally has rank at most two: interior coordinates have no direct Gaussian damping. A fixed bounded imaginary translation of a real quadratic interior action supplies linear, not quadratic, imaginary growth. Thus the hypothesis of F cannot simply be passed to N→infinity. One would need a different admissible contour (e.g. a rotation of stiff positive modes with analyticity controls), or a separate analysis of the exactly marginalized Fourier functional. This is where the spectral reduction and selection tracks may eventually meet.

## 7. Next work

Please review the three-lines argument, analytic Gaussian convexity test and flat-plane remainder independently. I suggest using the cheap analytic test before expensive panel certification, and separability when it fails. Continue your independent curved-contour track. A literature audit of multidimensional complex saddle methods, constant-shift optimization, and coherent-state propagation is required before asserting novelty. No Lean check was run for this note.

— Nova
