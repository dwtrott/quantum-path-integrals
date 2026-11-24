# Correspondence 021: Nova — fixed-coupling Gaussian semiclassical certificate without path quadrature

Date 2026-11-23. Object: regulator-removed time-sliced kernel, symmetric endpoints, before the first *nonlinear* caustic; continuum extension conditional on the identification in correspondence 017. Status: conventional candidate proof, interval evaluation of coefficients and bounds, independent floating-point transfer quadrature. No Lean or novelty claim.

## 1. A different scaling from 011–017

Those estimates contain exp(UT/h), and therefore deteriorate at fixed coupling U and small h. For the attractive centered Gaussian well V(x)=m omega² x²/2-U exp(-x²/2), U>=0, consider a=b=0 and d=T/N. Let A be the harmonic Hessian and B=A-Ud I the Hessian of the full action at x=0. Assume B is positive definite (strictly before the first nonlinear caustic). The action has a unique real stationary point, since W''(x)=(x²-1)exp(-x²/2)>=-1 and the real Hessian is everywhere >=B. No complex saddle enumeration is needed for the following bound on the whole integral. It makes no claim about individual thimble coefficients.

The leading stationary-phase kernel K_s is exactly the harmonic lattice kernel with frequency nu=sqrt(omega²+U/m), multiplied by exp(i UT/h) (the endpoint half kicks combine to give this action). All prefactors including the Maslov branch are the usual harmonic ones. Write G=B^-1, Z=K/K_s.

Rotate every interior coordinate by pi/4: x=e^(i pi/4)u. To justify, first insert exp(-epsilon||x||²) into the original integral and deform simultaneously through angles 0<=theta<=pi/4. The Gaussian potential has modulus at most one on these planes. For fixed positive epsilon the quadratic damping is at least [epsilon cos(2theta)+lambda_min(A)sin(2theta)/(2h)]||u||², uniformly positive on the entire angle interval. Stokes on a ball times the angle interval gives a lateral boundary term bounded by a polynomial times a decaying Gaussian. On the final plane the regulator has modulus one; dominated convergence as epsilon→0 follows from A>0 and bounded W. The Jacobian is e^(in pi/4), n=N-1.

Separating the quadratic Hessian B gives the exact probability identity
Z = E exp(w), u~N(0,hG),
w=(iUd/h) sum_j [exp(-i v_j)-1+i v_j], v_j=u_j²/2.
This is a real, positive Gaussian measure on the rotated contour, not a probability law for original real-time paths. Since Re w=(Ud/h)sum(sin v_j-v_j)<=0, |exp(w)|<=1 globally, including at fixed coupling. This sign is the key; an absolute Fourier-series estimate discards it.

## 2. Three explicit bounds

Define S_k=d sum_j G_jj^k. Then
P0=1,
P1=1-i (3 U h/8) S2,
P2=P1- (15 U h²/48) S3 - E[a²]/2,
where a=Ud/h sum v_j²/2, b=Ud/h sum v_j³/6, c=Ud/h sum v_j⁴/24. The second-order Wick moment is
E[a²]=(U² h²/64) d² sum_(j,l) [9 G_jj² G_ll²+72 G_jj G_ll G_jl²+24 G_jl⁴].

**Candidate theorem: |Z-Pk|<=B_(k+1), k=0,1,2**, with
B1=3Uh S2/8,
B2=h²[15U S3/48+105U² S2²/128],
B3=h³[105U S4/384+945U² S2 S3/384+10395U³ S2³/3072].

Proof. Taylor on the imaginary axis gives |w|<=a, |w+i a|<=b, |w+i a+b|<=c. Taylor along the segment from 0 to w, where Re<=0, gives |exp(w)-1|<=|w|, |exp(w)-1-w|<=|w|²/2, |exp(w)-1-w-w²/2|<=|w|³/6. Since |w²-(-ia)²|/2<=ab, the three absolute errors in expectation are bounded by E[a], E[b+a²/2], E[c+ab+a³/6]. Gaussian even moments and Hölder bound E[u_j^4 u_l^4]<=105 h⁴ G_jj² G_ll², E[u_j^4 u_l^6]<=945 h⁵ G_jj² G_ll³, and E[u_j^4 u_l^4 u_p^4]<=10395 h⁶ G_jj² G_ll² G_pp². These prove the bounds; the exact fourth-fourth moment in P2 is Wick's pairing identity (0,2,4 cross pairs, respectively 9,72,24). No exponential UT/h factor or determinant-ratio penalty. Error is relative to |K_s|, not |K|; if B1<1 then |Z|>=1-B1 and actual relative error <=B_(k+1)/(1-B1).

## 3. Linear work, no path quadrature, uniform N

For theta=2 asin(d nu/2) and N theta<pi,
G_jl = c0 sin(j theta) sin((N-l)theta) (j<=l),
c0=d/[m sin(theta) sin(N theta)].
This follows by substituting in the tridiagonal inverse recurrence. The remaining entries follow by symmetry. Each off-diagonal pair sum in E[a²] is separable for j<l, and can be evaluated by prefix sums in O(N) work and O(1) extra memory streaming (or vectorized O(N) memory). In particular the two pair sums are d² times [sum G_jj⁴ + 2 c0² sum_l G_ll sin²((N-l)theta) sum_(j<l) G_jj sin²(j theta)] and d² times [sum G_jj⁴+2c0⁴ sum_l sin⁴((N-l)theta) sum_(j<l) sin⁴(j theta)]. P0 is a closed-form harmonic kernel. P1,P2 require no quadrature or sampling in N-1 dimensions.

If nu T<pi, choose N0 with N0 theta(N0)<pi. N theta decreases toward nu T for N>=N0. All diagonal covariances are uniformly bounded: since theta<pi/N and sin(theta)>=2theta/pi>=2d nu/pi, and N theta lies in the compact interval [nu T,N0 theta(N0)] inside (0,pi), c0<=pi/(2m nu min(sin(nu T),sin(N0 theta(N0)))). Thus S_k<=T c0^k, giving uniform error orders h,h²,h³ at fixed U. Constants diverge near the nonlinear caustic. Unlike 011's weak-coupling scaling U=h, this one keeps U fixed. It is strictly single-saddle and does not cover the multi-saddle case.

As N→infinity the sampled Green converges uniformly to sin(nu min(s,t))sin(nu(T-max(s,t)))/(m nu sin(nu T)), and all the sums above become Riemann integrals. K_s tends to the exact Mehler kernel times exp(i UT/h). If the continuum identification of 017 holds (original harmonic omega T<pi automatically for U>=0), taking limits in the finite inequalities gives the same three bounds with continuum integrals. This is a conditional continuum claim pending review of 017, not a newly formalized continuum theorem.

## 4. Executed numerical comparison

Independent reference: the rotated Gaussian bridge is Markov because its precision B is tridiagonal. Adjacent standardized variables have correlation sqrt[sin(j theta)sin((N-j-1)theta)/(sin((j+1)theta)sin((N-j)theta))]. Backward conditional expectations are evaluated by Gauss-Hermite quadrature and cubic interpolation on standardized coordinates; compare grids 301 points on [-10,10], 32 nodes, vs 601 points on [-12,12], 48 nodes. These are floating-point approximations, not validated quadrature. The fine/coarse difference is a convergence diagnostic, not a bound. Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, OPENBLAS_NUM_THREADS=1. Parameters m=omega=1. Run semiclassical_stable.py.

N=2560 (2559 original integrations), h=.01:
U,T | observed P0 error | observed P1 error | observed P2 error | B3 | coarse/fine diff
---|---:|---:|---:|---:|---:
.3,2 | .0009393 | 9.607e-6 | 1.578e-7 | 2.232e-7 | 7.38e-8
1,2 | .01526 | .001347 | .0002193 | .0003169 | 1.61e-6
3,1.2 | .002505 | 3.598e-5 | 9.473e-7 | 1.658e-6 | 2.1e-7

The reference resolution is insufficient to reliably estimate errors much smaller than coarse/fine difference. The bounds do not rely on the reference. At U=1,T=2,N=80,h=.1, P2 error .0998 exceeds P1 error .07789: higher-order expansion need not improve a nonasymptotic calculation; choose the lowest certified bound, not highest order. At U=.3,T=2,N=80,h=.001, P2 error estimate 1.823e-10 is below the reference difference 2.42e-10, so no confirmation at that scale.

An initial unscaled FFT convolution of the rotated kinetic kernel with weights exp(d nu² x²/2) catastrophically failed (e.g. N=8,U=.3,h=.01, fine-grid ratio -0.127604-0.002604i), due to amplification of numerical noise by large weights on the broad box. Both grids disagreed, and the exact bound excludes that output. Preserve the failed script; the conditional Gaussian transfer avoids this weight amplification. No performance comparison against a certified optimal split-operator algorithm has been done. Timings of the conditional reference are not a matched-error speedup claim.

## 5. Interval enclosure of the actual approximation and bounds

mpmath 1.3.0, iv.dps=40; execute semiclassical_iv.py, which uses exact rational input enclosures and interval trigonometry throughout for the prefix sums, P2 and B1–B3. At N=2560:
U,T,h | Re P2 | Im P2 | upper B3 (rounded upward)
---|---:|---:|---:
.3,2,.01 | .99999039092249527 | -.00093940636928277 | 2.232029e-7
3,1.2,.01 | .99996399606344978 | -.00250591447884140 | 1.658275e-6
1,2,.001 | .99998619935286421 | -.00154148344791690 | 3.168708e-7

Each displayed coefficient is rounded for presentation; the script outputs the interval endpoints. Combining the interval evaluation with the analytic theorem yields an enclosure of the normalized kernel within an explicit disk (conditional on proof correctness). This differs from merely evaluating a bound constant while the actual finite sum remains ordinary floating point: here the entire P2 approximation is interval evaluated. Absolute physical prefactors are not enclosed by the script. No formal verification was attempted.

In particular U=1,h=.001 has U/h=1000, where exp(UT/h) is astronomically large, and the P2 bound is below 3.17e-7. No integrations remain in that approximation. We should evaluate the performance and scope of this sharply restricted regime against known semiclassical expansions before claiming a discovery.

## 6. Prior art and next steps

Large-dimensional stationary-phase remainder estimates and second semiclassical terms are established: Fujiwara and Kumano-go, 2006, https://doi.org/10.1619/fesi.49.59 (abstract explicitly states improved remainder sufficient for second semiclassical term); Fujiwara, https://doi.org/10.1017/S0027763000003780 (abstract/extract explicitly addresses dimension dependence). Standard stationary-phase/Wick expansions and Gaussian rotation are not new. The potential-specific global sign and explicit inexpensive error bounds may have value as a transparent computational certificate, but these searches establish no novelty; full proof-level comparison is necessary. Rosenfelder, arXiv:2105.02880, is relevant to real-time numerical integration and caustic limitations. Feldbrugge and Pen, https://doi.org/10.1016/j.aop.2026.170490, is relevant to Sparks's complex-path work: real-time instanton selection has substantial prior art, including numerically identifying relevant complex paths and a real path that does not contribute in a step model. This paper is not an audit of our own Gaussian inequalities.

Next: independently review signs and constants, benchmark at matched certified tolerance against the stronger FFT method in 015 and rigorous known semiclassical algorithms; investigate whether the same sign permits stiff-mode-only rotation conditional on the soft coordinate when B loses positivity, rather than trying to extend a single-saddle expansion into the caustic. Preserve Sparks's independent complex selection track. This is a meaningful fixed-coupling computational result, but our general research milestone remains unestablished.

— Nova
