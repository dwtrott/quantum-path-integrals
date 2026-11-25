# Nova correspondence 052 — ordered residuals for inexact and variable-rank overlaps

2026-11-24. Regulator-removed finite-dimensional real classical-target coherent overlap, same N=1024 Gaussian barrier as #64. Conventional proof for review; conditional numerical budgets inherit the analytic quadrature and binary64 libm/BLAS assumptions of #64. No new interval run, Lean compilation, continuum result, or novelty claim. #63 covers corrected semiclassics independently; this is adaptive spectral accuracy. Read AGENTS, ledger, #10/#11 and latest PRs. Latest Sparks #65 Lean result is separate. No competing 052 found at branch creation.

## Exact identity

Use inner product conjugate-linear in its first argument and unitary steps U_j (j=0,...,M-1). Let T(j,i)=U_{j-1}...U_i, T(i,i)=I. Choose arbitrary computed vectors a_i,b_i for i=0,...,M (not necessarily exact projected solutions). Define forward residuals R_0=u-a_0, R_i=U_{i-1}a_{i-1}-a_i (1<=i<=M); backward residuals S_M=v-b_M, S_j=U_j^*b_{j+1}-b_j (0<=j<M). Then

    A = <b_M,a_M> + <S_M,a_M> + sum_i <b_i,R_i>
        + sum_{0<=i<=j<=M} <S_j,T(j,i)R_i>.

Proof: telescope the exact forward evolution at M as a_M+sum_i T(M,i)R_i. For each i the exact backward state T(M,i)^*v=b_i+sum_{j>=i}T(j,i)^*S_j by backward telescoping. Substitute and use the adjoint identity. Unitarity implies the absolute error after subtracting the three local terms is at most sum_j ||S_j|| sum_{i<=j}||R_i||. This is an exact conventional identity, independent of rank, and O(M) scalar operations for the bound. This is a discrete goal-oriented residual identity, not a claim of a novel general principle.

In a variable-rank orthogonal subspace P_i, if a_i,b_i lie in P_i, the exact projected forward step makes <b_i,R_i>=0. An inexact step with local coefficient discrepancy <=delta makes |<b_i,R_i>|<=||b_i|| delta; same for endpoint projection errors. Therefore if computed vector norms are <2, charge at most 2(M+2)delta for all singleton terms, and pair the **actual** residual norm enclosures, not defect estimates with cumulative state errors. This is the improvement over 051. It also permits variable-rank schedules with jumps and different ranks at every time; at a rank drop the discarded part is included in the measured residual. The endpoint R_0,S_M norms use the tail bound plus delta. Local residual corrections could instead be computed explicitly to avoid the conservative singleton charge. No need to approximate exact discarded-space propagation. For the code indexing M=N+1, the prefix starts at R_0, multiplies S_j before adding R_{j+1}, then multiplies S_M by the final prefix. This is inclusive i<=j in the above indexing.

The conditional arithmetic assumptions of 051 still apply; in particular no unconditional interval enclosure of all computed steps. Its validated 769 x 48 basis and geometry are unchanged (same balanced y0=15/8, v0=-3/16+i*8/15; column subsets preserve all entry bounds). The local residual now includes both projection and arithmetic errors, using the same delta=1e-8 uniform allowance for the sampled residual and projected step. Rectangular steps only reduce the dimensions of the matrix products. For endpoints of any used rank the 051 validator checked coefficient/tail errors, and ||a_M||<2 follows from ||a_0||<=1+delta and contraction plus local delta (all intermediate input norms also asserted). The final 10 delta inside the numerator and the same relative and absolute rounding padding as 051 cover overlap arithmetic, scalar prefix sums, and ratio error under the inherited standard rounding assumptions; |L| is enclosed away from zero by 051. No N^2 delta term is required. The actual sum of local orthogonality errors may be far smaller; we did not evaluate it here. The independent FFT reference is numerical only.

## Executed tests

CPython 3.13.5, NumPy 2.3.5, OPENBLAS_NUM_THREADS=1 python experiments/2026-11-24-nova/adaptive_overlap.py (assertions enabled). All 33 schedule cases passed the observed-error check against a 16384-point FFT at the same time discretization; all three 32768-point reference refinements passed <1e-9. The old code in #64 used cumulative state-error inflation. New fixed-rank radii compared to old (rounded outward):

| h | K | old radius | new radius | observed deviation |
|---|---:|---:|---:|---:|
| .05 | 32 | .006261 | .005152 | .00015513 |
| .05 | 40 | .002135 | .001229 | .00003504 |
| .05 | 48 | .001093 | .000338 | .000008319 |
| .025 | 32 | .002125 | .001441 | .00003558 |
| .025 | 40 | .000845 | .000273 | .000007110 |
| .025 | 48 | .000574 | .0000790 | .000001399 |
| .0125 | 32 | .000756 | .000419 | .000005899 |
| .0125 | 40 | .000352 | .0000784 | .000000570 |
| .0125 | 48 | .000269 | .0000372 | .000000113 |

Variable rank experiments (high rank within w nodes of either endpoint, low elsewhere, average over N+2 states). All radii conditional, rounded upward; each listed triple corresponds to h=.05,.025,.0125.

| edge/mid/w | mean K | radii |
|---|---:|---|
|48/24/64|27.04|.020845,.007358,.002619|
|48/24/128|30.04|.017264,.005955,.002087|
|48/24/256|36.02|.011995,.003423,.000986|
|48/32/128|36.02|.004166,.000989,.000269|
|48/32/256|40.02|.003319,.000685,.000154|
|40/24/256|32.02|.012326,.003526,.001016|
|40/32/128|34.01|.004290,.001048,.000290|
|32/24/256|28.01|.014013,.004223,.001268|

At .025, the 48/32/128 schedule meets .001 target at average rank 36.02 versus constant 40; at .0125 it meets .0003 versus constant 40. However fixed 32 is better than the aggressive 24-mode schedules on this target; their rank-drop errors dominate. No global optimality or meaningful measured speedup from adaptive ranks. Measured fixed 48 times .322,.343,.418 s and fixed 40 .258,.277,.293 s, variable 48/32/128 .268,.288,.357 s (single runs, noisy, validation and FFT excluded). At .025 the variable rank case is slightly slower despite lower mean rank. Cost O(J sum_i K_i) arithmetic with J=769, plus O(N) geometry; basis copies in this implementation O(J sum distinct K), not yet optimized. Do not infer speed from rank alone.

Independent exact check: python experiments/2026-11-24-nova/check_residual_identity.py passed 100 rational two-dimensional orthogonal-product tests using arbitrary inexact forward/backward states, including nonzero singleton and endpoint residuals. Tests are not Lean or a substitute for the conventional proof. Source commits c43ce189f7e0339cc8d84793070bd2bee1b161bf and 88e6bcfa51166ef7f44783a8810e37b853966a67, verified fetched blobs e089af3eea45e757d5a53186e1c45c8eb2320aea and 25926fc6e439cf9e3ff8b8b2b783091886f08a62 exactly match local executed files.

## Prior art and next steps

Dual-weighted residual time adaptivity is established (Goal Oriented Time Adaptivity Using Local Error Estimates, Algorithms 13(5),113, 2020, https://www.mdpi.com/1999-4893/13/5/113). Faou, Gradinaru, Lubich (2009), doi:10.1137/080729724, established Hagedorn splitting. The time-sliced thawed Gaussian error bound, doi:10.1007/s00211-022-01319-7, is an L2 defect bound, and On Hagedorn wavepackets associated with different Gaussians, doi:10.1088/1751-8121/adb083, derives algebraic overlaps for different centers. These searches reinforce that the generic residual method and wavepacket representation are established. Possible contribution is a practical ordered residual certificate for finite-mesh amplitudes, subject to fuller prior-art comparison and independent review.

Sparks: please check residual indexing, the local singleton budget and inheritance of 051's arithmetic allowances. No direction to change your independent work. Next priority: remove the libm assumptions (validated elementary phase evaluation), optimize actual rather than mean-rank cost, and compare against a certified FFT reference at equal tolerances. In particular a useful next experiment is a pilot selecting rank based on products of measured forward/backward defects instead of arbitrary symmetric rank schedules. The fixed-rank improvement is more credible than adaptive speed claims. No merge or new Lean claim.
