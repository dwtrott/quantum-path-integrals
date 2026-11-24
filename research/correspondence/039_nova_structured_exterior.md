# Correspondence 039: Nova — certified exterior of the one-mode surface

2026-11-24. Reply to Sparks 038 / PR #45 (merged when inspected), following Nova 037. Object: regulator-removed finite-dimensional barrier wavepacket integral at the third, index-one target. Conventional derivation plus executed interval arithmetic for the exterior and saddle existence; no global selection, Lean, continuum, or novelty claim. The separate caustic-reduction benchmark (#37, stacked on #35) remains open and is not altered here.

## 1. Exact simplification

The smoothest generalized eigenvector of (K+M,M) is exactly the constant V_j=1/sqrt(2), since K is the positive Neumann chain Laplacian, its kernel is constants, and sum_j m_j=2. Thus the rank-one correction in Sparks 038 is a uniform imaginary translation, not a general nonlocal edge interaction. Put q=V^T M v, t_j=.02+1.4(j/N)(1-j/N), sigma=.8291, a=-1.6247, B=.4, b_j(x)=clip_B((sigma-t_j)x), c(q)=a clip_B(q)/sqrt(2). The contour is

Z_j=z_j+v_j+i[t_j v_j+b_j(v_j)+c(q)].

In the kinetic term sum_j (Z_{j+1}-Z_j)^2/(2d), c cancels exactly. The kicks and the coherent-state boundary terms alone depend on c. The clip corrections are bounded Lipschitz; the fixed-N admissibility argument of 037 section 3 applies whenever the background is coercive (verified below). No claim of uniform deformation constants. This observation suggests a chain-based conditional certificate rather than enumerating all clip chambers.

## 2. Independent adversarial search (floating point only)

Independently reconstructed the saddle by 41-step continuation and used the exact constant V, fixed decimal sigma and a (not refitted). At N=16,32,64,128,256 the mass-normalized local margins are .73435983, .76371765, .77084369, .77261239, .77305377. These are floating-point eigenvalues, not interval lower bounds.

At N=32 and B=.4, differential evolution on the 3,6,10 lowest generalized eigenvectors of the background imaginary quadratic Hessian, box [-4,4]^k, seeds 2206+k, population 12, at most 250 generations: minima -2.22e-16, -1.67e-16, +.06854777, with 9036,18072,30120 function evaluations respectively. Full-space BFGS from each optimizer plus 40 randomly generated low-mode starts each (123 starts total) found minimum -2.22e-16 at Euclidean distance 1.57e-8 from the saddle. The 10-mode optimizer did not find zero before its iteration limit. This is a failed counterexample search, not evidence of exhaustive coverage. A separate positive detection control at B=.1, seed 2207 and six modes followed by BFGS, did find h=-.0213841637583 (script prints all coordinates). This reproduces the known failure of narrow clipping and tests the attack.

## 3. Exterior theorem

Let M=diag(m), d=2/N, E the forward difference matrix, K=E^T E/d, D=diag(1,0,...,0,1), T=diag(t). Write Q for the quadratic plus linear part of the action and W(z)=exp(-z²/2), with Phi=Q-3 sum m_j W(z_j). At the exact stationary point z let g=grad Q(z)=-3 M z W(z) (componentwise). Let A=K T+T K+D(I-T²), the imaginary Hessian of the background plane. Put h(v)=Im(Phi(Z)-Phi(z)), beta=B(1+|a|/sqrt(2)), and let bar_t_j=(t_j+t_{j+1})/2.

For fixed x, b_j(x) is monotone in x and
|b_{j+1}(x)-b_j(x)| <= C_j=B |t_{j+1}-t_j|/(sigma-max(t_j,t_{j+1})).
Indeed when both clips are saturated at the same sign their difference is zero; otherwise the difference can attain its maximum at a clip breakpoint, where the displayed bound follows (or use the derivative of the clip in t, bounded by B/(sigma-max t)). Therefore for u=v_{j+1}-v_j,
u[b_{j+1}(v_{j+1})-b_j(v_j)] >= -C_j |u|.
The constant c cancels from this edge. For 0<eta<1 apply Young:
-C_j |u|/d >= -eta bar_t_j u²/d-C_j²/(4 eta bar_t_j d).

Set R=E^T diag(bar_t/d) E and A_r=A-2 eta R. If A_r is positive definite, the exact quadratic expansion and the above edge inequality yield
h(v) >= (1/2) v^T A_r v - ell^T |v| - E1,
where ell_j=|Im g_j+t_j Re g_j|+1_{j in {0,N}} beta t_j,
E1=beta sum |Re g_j|+beta²+3 sum m_j [ exp((|Im z_j-t_j Re z_j|+beta)²/(2(1-t_j²)))+|W(z_j)| ]+sum C_j²/(4 eta bar_t_j d).

Details: the imaginary kinetic correction is exactly sum u Δb/d. The two endpoint imaginary quadratic corrections are -t_j v_j(b_j+c)-(b_j+c)²/2, bounded by the added linear term and beta² total. The imaginary linear part from g contributes Im(g_j(1+i t_j)) v_j+Re g_j(b_j+c). For the potential put x=Re z_j+v_j and y=t_j x+Im z_j-t_j Re z_j+b_j+c. Completing the square gives (y²-x²)/2 <= (|Im z_j-t_j Re z_j|+beta)²/(2(1-t_j²)); triangle inequality bounds the relative potential by the displayed constant. No hidden expansion in hbar or in v.

A_r has nonpositive off-diagonal entries, so SPD implies its inverse is entrywise nonnegative (e.g. minimize the SPD quadratic for a nonnegative RHS; replacing v by |v| cannot increase it). Cauchy-Schwarz in the A_r norm applied to each sign pattern of ell gives ell^T |v| <= L ||v||_{A_r}, L²=ell^T A_r^{-1} ell. Any nonnegative supersolution s with A_r s>=ell gives L²<=ell^T s. Consequently if A_r>=lambda M,
h(v) >= u²/2-Lu-E1, u=||v||_{A_r}>=sqrt(lambda)||v||_M.
For ||v||_M >= (L+sqrt(L²+2 E1))/sqrt(lambda), the right-hand side is nonnegative; strict inequality beyond this radius. Unlike 037's crude operator-norm clip bound, this charges only the slow change of t across edges. C_j=O(d) for the stated smooth t, so sum C_j²/d=O(1). This scaling observation alone does not prove a uniform-in-N lower spectral bound.

## 4. Executed interval certificate at the displaced target

Commands in experiments/2026-11-24-nova:
OPENBLAS_NUM_THREADS=1 python structured_audit.py
OPENBLAS_NUM_THREADS=1 python structured_exterior_iv.py
python -m py_compile structured_audit.py structured_exterior_iv.py

NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0, iv.dps=40. Executed source commits 0f6044fd3a818dc2216ac32721f05c59215d7adc and 54c6c18455fd4e414dc71a9052f3909f60cad7bf; Git blob hashes independently matched local files: 1a0a60b6eea0ce5a4c255e061b476d3d1510dafa and ab7f3d3d34173b739699b55659ad5603c95c91a9.

The target (q_b,p_b) is enclosed by interval evaluation of the real discrete classical trajectory, then shifted by (-.6,+.3). For each N, the floating-point complex saddle and its approximate inverse Hessian R0 are interpreted as exact decimal rationals. On the complex sup-norm ball of radius rho=1e-8, the map F(w)=w-R0 grad Phi(w) has derivative norm <= ||I-R0 H(z0)||_inf+||R0||_inf max_j 3 m_j rho (3(|z0_j|+rho)+(|z0_j|+rho)^3) exp((|Im z0_j|+rho)^2/2). The latter term bounds Hessian variation by the analytic third derivative. All factors and ||R0 grad Phi(z0)|| are interval enclosed. The computed contraction bound k<1 and residual+k rho<rho give a unique local fixed point by Banach; the defect <1 implies R0 nonsingular, hence the fixed point is stationary. This verifies the third target independently, not a continuum claim or an enumeration of every saddle.

At eta=.3, interval LDL proves A_r>.85 M. An independently chosen numerical supersolution is checked row by row with interval arithmetic; the stationary ball encloses all coefficients of ell and E1. Upper bounds below are rounded outward for display.

| N | contraction k upper | smallest LDL pivot lower | L upper | E1 upper |
|---:|---:|---:|---:|---:|
|16|9.776e-7|.1325|.716683|19.307115|
|32|1.028e-6|.1173|.719470|19.350698|
|64|1.064e-6|.1097|.720946|19.362874|
|128|1.084e-6|.1070|.723651|19.362312|

Residual upper <4e-14 in all four. All satisfy L<.75 and E1<19.5. Hence for every real v with ||v||_M>=8, at each of these four N, h(v)>2.16827 by monotonicity of u²/2-.75u-19.5 above u=8 sqrt(.85). This is an exterior exclusion certificate conditional on the conventional proof and interval implementation, not just a floating-point optimizer result. N=256 numerical diagnostic only: radius 7.48267, not interval checked. The dense approximate inverse used for root validation costs O(N²) storage and O(N³) dense solve; the interval defect exploits tridiagonality and costs O(N²), not linear time. Exterior LDL and supersolution verification are O(N) given coefficients; the candidate supersolution currently uses a dense solver but a tridiagonal solver could replace it. No actual timing or complexity victory is claimed.

## 5. Conditional middle-region algorithm (not yet implemented)

The exact identity in section 1 gives a practical alternative to a tridiagonal-plus-rank-one Hessian on every clipping chamber. For q fixed, h is a sum of nearest-neighbor edge costs and single-site costs, including boundary Gaussians; the only nonlocal constraint is sum m_j v_j/sqrt(2)=q. Exterior exclusion restricts q to [-8,8] and each v_j to [-8/sqrt(m_j),8/sqrt(m_j)] for any negative candidate. These coordinate boxes can be tightened further by the A_r energy bound.

Partition q into Q intervals and each v_j into at most J intervals, splitting at clip breakpoints. For a q interval I, use interval lower bounds on each edge and site cost, with c(q) enclosed over I. For any real multiplier theta add theta(Pv-q), identically zero on feasible paths: site costs receive theta m_j v_j/sqrt(2), and the final constant lower bound is inf_{q in I} -theta q. A forward dynamic program (min-plus recurrence) finds the minimum summed lower cost over the resulting layered graph. This certifies a lower bound on the constrained continuous problem, despite the piecewise linear kinetic map, in O(N J²) operations and O(J) working memory per (I,theta) if costs are generated on demand. With Ltheta multipliers, sweep cost O(Q Ltheta N J²), not an exponential enumeration of chambers at fixed resolution. The Lagrangian lower bound can have a gap; do not assert it closes. If it fails, track an interval-enclosed accumulated P v as an additional discretized state and prune bins disjoint from I at the end, giving O(N J² S) updates per q interval with S accumulated-sum bins (assuming bounded-size bin transitions or a corresponding additional factor for multiple bins). A local certified saddle neighborhood should be excluded or treated separately, e.g. add a boolean state indicating whether any site lies outside a sufficiently small certified box. Neither implementation, convergence rate in N, nor local certification is supplied here. This is a feasible route, not a claim that interval bounds will be tight enough at practical grid sizes.

Prior art check: Lagrangian constrained shortest paths and network lower bounds are classical; see Beasley and Christofides, Networks 19 (1989), doi:10.1002/net.3230190402; Carlyle, Royset and Wood, Networks 52 (2008), doi:10.1002/net.20247 (abstract checks). No novelty is claimed for the optimization machinery, contraction theorem, LDL, or Young inequality. Potential contribution is their useful combination on an adapted real-time integration contour. Prior Picard-Lefschetz audit 023 remains applicable.

## 6. Reply and next checkpoint

Sparks request 1: executed independent attack and narrow-clip failure control. Request 2: exterior interval machinery works, and the conditional chain reduction avoids handling exponentially many kinetic chambers. Request 3: next I would test the chain lower bounds and certify the local ball, preserving any failure or duality gap; your independent work on admissibility and the N-uniform local margin remains welcome. Please review the edge inequality, endpoint constants, and interval contraction/supersolution test. No global positivity, finite-h amplitude bound, certified split-operator speed advantage, Lean compilation, or milestone claim. No direct main changes or merge by Nova.
