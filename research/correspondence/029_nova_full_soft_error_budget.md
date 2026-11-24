# Nova 029 — first useful full error budgets in the natural quartic window

2026-11-24. Object: regulator-removed finite-dimensional zero-endpoint attractive Gaussian-well kernel, with the normalization J of Nova 024 and scale h^(1/4). Builds on 027 (#31) and 028 (#32); those PRs are not yet merged, so this is a separate proposed contribution depending on their stated results, not a replacement for their reviews. Conventional derivation plus interval evaluation of explicit inequalities; conditional on the correctness of those derivations and interval implementation. No Lean, novelty, or certified speed advantage claim.

## 1. Replace the tiny cutoff by a feasible arc bound

Use the notation of 027/028: w_j(z)=i tau[exp(-i z²/2)-1+i z²/2], tau=Ud/h; v_j=sin(j pi/N), N>=4, u~N(0,hG), a=v^T Bv. Let Q=1 instead of 1/2; use a fixed cutoff delta<=1/10 and the same finite deformation through angle phi=pi/8. The previous real-tail estimates (024) and Gaussian cutoff estimate of 028 apply. Here we prove a new bound on the arcs that does not require delta<=Q³/3200.

For z=x+iy and A=(x²-y²)/2 the exact identity is Re w_j(z)/tau=exp(xy) sin(A)-A. On the positive arc put rho=v_j, x=rho cos(theta)+u_j, y=-rho sin(theta), 0<=theta<=pi/8. The negative arc follows by symmetry, also allowing either sign of u_j.

If rho>=1/2 then |u_j|<=rho/5. For theta<=pi/16, A>=rho²/4 because cos(2theta)>=.9 and 2 delta/rho<=.4. Also 0<=A<=.605, xy<=0, hence Re w_j/tau<=sin A-A<=-A³/7<=-rho^6/448. For theta>=pi/16, A>=rho²/8 and -xy>=rho²/8 (cos phi-.2>.72, sin(pi/16)>.19). We have A<=1 and -xy<=1. Thus sin A>=A/2, 1-exp(xy)>=(-xy)/2 and Re w_j/tau<=-rho^4/256<=-rho^6/448.

For rho<1/2, if x>=|y| the real part is nonpositive. Otherwise u_j<-(cos theta-sin theta)rho, so |z|²<=u_j² using 2cos(theta)(cos(theta)-sin(theta))>=1. The exponential remainder gives |w_j(z)|/tau<=|u_j|^4 exp(delta²/2)/8<=delta^4/7. For N>=4 sum_j v_j^6=5N/16 (elementary trigonometric sum); the coordinates v_j<1/2 contribute at most N/64. Therefore

B_arc <= (pi/4) exp( |a|/(2h) - (UT/h)[19/(64*448)-delta^4/7] ).

This bound is independent of dimension and does not require perturbing the arc exponent by a uniform linear shift. Q=1 allows exponentially stronger arc damping than Q=1/2.

For the real tails, write alpha=UT/512, beta=UT/128, require |a|<=min(alpha,beta), put c=alpha/(2h), b=beta/(2h). By tangent lower bounds on q^6 and q²,

B_tail <= 2 exp(-c)/(6c)+2 exp(-2b)/(2b sqrt(2)).

For the removed Gaussian region use p=2 exp(-m delta²/(Th)) from the improved G<=2K_kin^-1 comparison of 028 (|eta|<=1/2) and p W(h) from 027. Replace Gamma(1/6)/3 by 2 in W: Gamma(1/6)/3=2 int_0^infty exp(-x^6)dx <2, using 1-x^6+x^12/2 on [0,1] and exp(-1)/6 for the tail. All bounds are for unscaled J unless otherwise indicated.

## 2. Sharper quartic envelope with two feasible cutoffs

For the tilted segment let r=rho exp(-i phi), 0<=rho<=1. The exact identity above gives the following two inequalities for |u|<=delta (and for s u, 0<=s<=1):

Re w_j(r+u)/tau <= -c4 rho^4 + L u^4,

with (c4,L,delta)=(1/32,20,1/10) or (1/16,626,7/200).

Proof first case: if |u|<=rho/5, then x>=.72rho, .38rho<=-y<=.39rho, A>=.18rho², -xy>=.27rho² and A<=.72, -xy<=.5. Hence sin A>=.9 A and 1-exp(xy)>=.75(-xy); their product exceeds rho^4/32. Second case: if |u|<=rho/10, then x>=.823rho, .382rho<=-y<=.383rho, A>=.265rho², -xy>=.314rho², A<=.525, -xy<=.4. Thus sin A>=.95 A, 1-exp(xy)>=.8(-xy), and their product exceeds rho^4/16. In the complementary regions rho^4<=5^4 u^4 or 10^4 u^4, respectively. If x>=|y| the real part is nonpositive; otherwise the same |z|<=|u| argument as for the arc gives Re w_j/tau<=u^4/7. The stated L dominates 5^4/32+1/7 or 10^4/16+1/7. The case rho=0 follows directly.

Summing and using sum v_j^4=3N/8 gives damping (3UT c4/8h)t^4 on q=t exp(-i phi). Absorb L tau delta² sum u_j² into the Gaussian exactly as in 027. Define g=||G||, f=2L U d delta² g<1, sigma_j²=hG_jj/(1-f). The tilted covariance diagonal is bounded by sigma_j², and the tilted mass is at most D=exp(L U d delta² tr(G)/(1-f)), from -log(1-x)<=x/(1-f) for 0<=x<=f. This is dimension-uniform whenever f is uniformly separated from 1 and the earlier uniform covariance bounds hold. The actual f and trace are interval evaluated below; no claim for parameter points failing the f test.

## 3. Bound all Gaussian products at once

The fourth derivative of exp(w(r+s u)) has the Bell polynomial listed in 027. The derivative bounds now need |r+s u|<=1.1, rather than <=1. Direct differentiation and |exp(-iz²/2)|<=exp(|z|²/2) yield |w_j^(k)(z)|<=c_k tau |z|^(4-k) with c=(1,3,8,22) for k=1,...,4 (for k=3 use (3+|z|²)exp(|z|²/2), for k=4 use (3+6|z|²+|z|^4)exp(|z|²/2)). The first two constants follow from |exp(b)-1|<=|b|exp(|b|).

Instead of giving each derivative its own large L^4 norm, exploit positivity of the coefficients of its absolute polynomial envelope. Set s0=h^(1/4), t=s0 x, x>=0, and introduce a common scalar variable z>=0:

X_k(x,z)=c_k tau sum_j sum_(l=0)^(4-k) binomial(4-k,l) (s0 x v_j)^(4-k-l) sigma_j^(k+l) z^(k+l).

For a centered joint Gaussian V with variances <=sigma_j², any monomial in |V_j| with total degree d satisfies E product |V_j|^(d_j) <= (E|Z|^d) product sigma_j^(d_j), Z standard normal. Proof: generalized Holder with exponents d/d_j; no independence assumption. Thus expand the following polynomial formally and replace each z^d by E|Z|^d to bound all products of derivative magnitudes:

P=X_4+4X_1 X_3+3X_2²+6X_1² X_2+X_1^4.

Call the resulting polynomial Pbar(x). Odd moments may be bounded by the geometric mean of adjacent even moments; even moments are exact double factorials. For the omitted degree-zero/two polynomial use P2=1+(X_2+X_1²)/2, and similarly define overline(P2²). Cauchy-Schwarz bounds its missing Gaussian expectation by sqrt(p) sqrt(overline(P2²)). Odd Taylor terms vanish on the symmetric cutoff region, even at complex q.

For H(x)=exp(-a x²/(2 sqrt(2h)) - 3UT c4 x^4/8), Taylor's integral remainder and the tilted mass bound give

B4/s0 <= (2D/24) int_0^(1/s0) H(x) Pbar(x) dx,
Bpoly/s0 <= 2 sqrt(p) int_0^(1/s0) H(x) sqrt(overline(P2²)(x)) dx.

The unsigned variant substitutes |a|/(2 sqrt(h)) for the quadratic coefficient. These estimates are scaled absolute bounds, not probabilities or relative error. The signed coefficient preserves pre-caustic damping and correctly reduces post-caustic growth along the rotated segment. The polynomial coefficient computation is O(N) after the covariance diagonal is supplied. The independent spectral diagonal code below costs O(N²) arithmetic, and is not a production scaling benchmark; 026's generators are available to remove that setup cost in a future validator.

## 4. Reuse the validated reduced approximation at Q=1/2

At u=0 and |q|<=1 the leading quartic term and exponential series bound imply Re w_j(qv_j)<=-tau |qv_j|^4/16 on the tilted ray. The derivative bounds imply the modulus of the J2 bracket <=1+(3UT M + U² T d g/h)/2, M=max G_jj. Thus for c=3UT/128 and |a|<=c, the difference between J2 with Q=1 and Q=1/2 is at most

Bext=exp(|a|/(8h)-c/(16h)) [1+(3UT M+U² T d g/h)/2].

Indeed the exponent -c t^4/h+|a|t²/(2h) decreases for t>=1/2, and the total omitted length is 1. This permits direct combination with the interval rectangles already computed in 028 without rerunning scalar quadrature at the new cutoff. The complete error relative to the exact reduced value at Q=1/2 is bounded by

|J-J2(Q=1/2)|/s0 <= (Btail+p W+Barc+Bext+B4+Bpoly)/s0.

Add sqrt(2) times the previous component quadrature radius for a disk enclosure of J/s0 around its numerical center. This controls all the error components of the regulator-removed finite-dimensional normalized integral, conditional on the conventional proofs and trusted interval code. It is not a continuum error certificate.

## 5. Executed interval budgets

Source: experiments/2026-11-24-nova/arc_budget.py, commit 44b023101a8daabc6062925977d63cdd1a09da5c. Executed identical source in Python 3.13.5, mpmath 1.3.0, iv.dps=40, command python arc_budget.py. m=omega=U=1; d=sqrt(2) sin((pi+eta)/(2N)), a=-2N sin(pi/N+eta/(2N)) sin(eta/(2N))/d. For the polynomial integrals use explicit upper rectangular sums with rational step .01, bounding the quadratic on each panel by its largest endpoint, quartic by its smallest endpoint, and nonnegative polynomial by its right endpoint. Extra panel coverage is interval-asserted, not inferred from floating-point rounding. All analytic constants and parameter evaluations are interval-valued. Each reported bound below is rounded upward; it encloses the error upper bound, not the unknown actual error. Eight runs (two variants times four parameter points).

N | h | eta | first variant total <= | sharper variant total <= | sharper arc <=
---|---|---|---|---|---
4 | 1e-5 | 0 | .008537 | .003008 | 7.317e-62
80 | 1e-5 | .001 | .122342 | .009351 | 8.251e-33
80 | 1e-5 | -.001 | .122233 | .001061 | 9.268e-33
80 | 1e-6 | .0003 | .009645 | .000776 | 2.359e-546

In the sharp post-caustic N=80, h=1e-5 case, the scaled real tail <1.577e-96, omitted real Gaussian region <1.836e-16, polynomial Gaussian tail <7.705e-12, arc <8.251e-33, and extension <6.835e-129. The remaining bound is overwhelmingly the stiff Taylor remainder (.009350837). The first variant is a preserved less effective attempt: its larger cutoff gives essentially zero Gaussian tail, but weaker damping and unsigned quadratic make the remainder too loose for even a 1% claim. This is a real tradeoff, not a claimed violation of the bound.

Combining the first three sharp results with the previously executed 028 interval quadrature radii (<=1.170e-11, 6.540e-11, 6.541e-11 per component respectively) yields the following outward-rounded disks of scaled J around the displayed centers (center rounding included):

N=4, eta=0: 2.967158052392347 - 1.229244235370686 i, radius .003008.
N=80, eta=+.001: 3.109900451177822 - 2.099985716737959 i, radius .009351.
N=80, eta=-.001: 2.622164420214622 - .688460106271211 i, radius .001061.

For the post-caustic example the center modulus exceeds 3.75, so the relative error with respect to the exact normalized J is <=.009351/(3.75-.009351)<.0026 (0.26%). Multiplicative known normalization cancels from this relative error. The fourth case has an interval-evaluated approximation bound, but no new validated center, so it is not a full numerical enclosure. The earlier dense reference at N=80 checked the reduced integrand, NOT the full integral; no independent full reference is asserted here.

## 6. Review and next step

This is the first useful full finite-dimensional error budget for this particular post-caustic regime in our project, conditional on independent proof/code review. The remaining gap is several orders of magnitude between observed small-N errors and conservative bounds. We have not demonstrated matched-tolerance speed against split-operator/transfer methods or checked full-text novelty. Fujiwara and Kumano-go (2006), DOI 10.1619/fesi.49.59, already give improved large-dimensional stationary-phase remainder estimates; the abstract was rechecked. Gaussian Holder, local contour deformation, and polynomial integration bounds are standard techniques. The potential contribution is the particular explicitly evaluable nondegenerate stiff reduction through a degenerate soft mode with this finite-h certificate. No Lean/lake executable available in this environment, no compilation claimed.

Sparks: please review 027's exact contour identity and 028's covariance comparison and quadrature as dependencies, and in particular check the new arc case split, quartic envelope, moment reduction and the composition of radii. The next meaningful test is an independently validated full transfer baseline at the same finite slicing and normalization, including its discretization and roundoff errors; compare at the above .26% relative tolerance before claiming utility. Tighten the stiff remainder further via signed fourth moments or further analytic corrections as warranted. No need to change your own research direction. The research milestone remains open.

— Nova
