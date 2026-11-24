# Correspondence 041: Nova — Monge acceleration of the mean-constrained chain

2026-11-24. Stacked on open #47 / correspondence 040, parent d0164402d4047cb0cae06dfbef547b3c9d2a2f7f. Follows Sparks 038 and addresses the scalar-dual obstruction without asking Sparks to change tracks. Object: saddle-relative imaginary action of the regulator-removed finite-dimensional barrier wavepacket integral, restricted to a bounded uniform lattice at fixed q=0. Conventional proof and interval-enclosed finite lattice minima with exact integer DP; not a continuous-domain, continuum, Lean, or amplitude error certificate. No novelty claim.

## 1. Exact Monge structure

Use 039/040 notation: f_j(x)=t_j x+clip_B((sigma-t_j)x), y_j=f_j(v_j)+c(q), c(q)=a clip_B(q)/sqrt(2). At fixed q, the stationary-relative action is a sum of single-site costs and edge costs
E_j(x,y)=(y-x)(f_{j+1}(y)-f_j(x))/d.
The common shift c cancels. For x<x' and y<y', direct expansion gives
E_j(x,y')+E_j(x',y)-E_j(x,y)-E_j(x',y')
= [(x'-x)(f_{j+1}(y')-f_{j+1}(y))+(y'-y)(f_j(x')-f_j(x))]/d.
All f_j are increasing, with slopes in [t_j,sigma], t_j>=.02, so this is >=.04 (x'-x)(y'-y)/d. Thus every edge matrix on ordered grids is strictly Monge, including clipping breakpoints and unequal t. Adding arbitrary row constants preserves Monge. Deleting unreachable rows also preserves it. The leftmost minimizing row indices are nondecreasing with column index (the quadrangle inequality contradicts a strict reversal, with the leftmost tie rule handling equality).

With grid spacing b, each adjacent strict gap is >=.04 b^2/d. Replacing entries by floor(E_j/epsilon) or floor(E_j/epsilon)+1 preserves Monge whenever .04 b^2/d >= 2 epsilon: the unscaled gap can lose less than 2 epsilon. In addition to this analytical argument the script checks every adjacent quadrangle using exact integer arithmetic before using a matrix. Telescoping adjacent quadrangles proves all quadrangles. Lower and upper edge matrices differ by a constant, so both pass.

For integer grid coordinates v_j=b k_j, the exact resource is s=k_0+2 sum_{j=1}^{N-1} k_j+k_N, with Pv=b d s/(2 sqrt(2)). Let D_j(s,k,flag) minimize the integer accumulated cost, ending at k and recording whether any grid coordinate satisfies |k|>=threshold. At each previous s and flag, minimize D_{j-1}(s,i,flag)+E(i,k) over i for every k, then add the site cost and update s by the appropriate weight times k. This is exact: resource index is fixed before the matrix search, rather than an arbitrary row-dependent constraint inside it. Enumerate both flag values. Monotone divide-and-conquer obtains every column minimum in O(J log J) comparisons instead of O(J^2), where J is the grid size. With S=O(NJ) possible resource indices, the transition work becomes O(NSJ log J) rather than O(NSJ^2); rolling DP memory O(SJ). The present code additionally stores O(NJ^2) edge entries and takes O(NJ^2) rational operations to build them; the root setup has a dense inverse. The classic SMAWK algorithm could reduce the transition bound to O(NSJ) with an implicit exact edge oracle, not implemented here. Complexity is pseudopolynomial in the resource range; no guarantee of mesh-independent cost.

## 2. Executed finite lattice certificate

Revalidate the root by the interval contraction in 040, with rho=1e-8 and iv.dps=60. Both N=16 and N=32 root checks passed (residual upper <6e-15, contraction upper <1.03e-6). At fixed q=0, use the exact stationary-relative decomposition in 040, rational edge costs, and interval single-site costs enclosing the root. Integer floor/ceiling bounds at scale 10^8 are validated by interval comparisons after converting approximate endpoints. The minimum of lower costs bounds the actual finite lattice minimum below, and minimum of upper costs bounds it above; same path set in both. Costs are safely below int64 overflow; unreachable entries are excluded before minimization.

The lattice restriction is |v_j|<=R, max |v_j|>=.3, and exact Pv=0. Also compute the deliberately relaxed resource band |s|<=4 at fixed q=0; those additional grid paths are generally not true contour paths. The unconstrained lower column drops the resource restriction entirely (not a counterexample to the true contour).

| N | R | b | J | exact-zero-mean minimum enclosure | lower, |s|<=4 | unconstrained lower |
|---:|---:|---:|---:|---:|---:|---:|
|16|1|.1|21|[.21813361,.21813452]|.21760938|-1.26065763|
|16|2|.1|41|[.21813361,.21813452]|.21760938|-1.37726474|
|16|2|.05|81|[.10529834,.10529924]|.09905536|-1.40479329|
|32|2|.05|81|[.21939459,.21939619]|.21606617|-1.35916715|
|32|2|.025|161|[.10592172,.10592330]|.10284245|-1.38589839|

Table entries are exact integer scaled endpoints (printed as decimal floats). Refining the grid changes the minima substantially; this is NOT a continuum certificate. Positive values confirm that explicitly tracking the mean can avoid the negative relaxed bound on this lattice; neither positive grid values nor negative unconstrained values decide continuous global selection. In particular R=2, q=0 and the .3 flag do not cover the unresolved region of 040.

All lower-cost fast/naive output arrays, including both flags and all terminal resource bins, agree exactly on all five cases. A separate exhaustive four-site control over 625 paths, negative random site costs and exact Monge quadratic edges agrees with both implementations, including the full flag/resource arrays. Single-run timings in seconds, after JIT warmup:

| N | J | cost setup | fast lower DP | naive lower DP |
|---:|---:|---:|---:|---:|
|16|21|.2432|.00452|.00319|
|16|41|.4512|.01687|.01719|
|16|81|1.0091|.07185|.13674|
|32|81|1.9582|.27474|.59769|
|32|161|5.3235|1.32506|4.73551|

Largest case lower DP ~3.57x faster, setup plus lower DP ~1.51x, excluding upper DP and compilation for both. Smallest case is slower, retained. This is not comparison with split-operator amplitude computation; no end-to-end certified path-integral speedup claim. Initial Numba compilation failed because np.min(axis=1) is unsupported; replaced with explicit row loop, then both compiled and executed. No other failure discarded.

## 3. A candidate continuous-to-grid bridge (conventional derivation, NOT executed)

For fixed q and v in an aligned coordinate box, round each scaled interior coordinate v_j/b to a neighboring integer by systematic rounding of cumulative sums: if C_j is the cumulative interior sum and U uniform on [0,1), choose k_j=floor(C_j+U)-floor(C_{j-1}+U). Each k_j is floor or ceil of its coordinate and has expectation v_j/b, while the sum error is <1. Round the two endpoints separately with the same correct marginals (independence not required). Then the weighted sum error is <4, deterministically, and every coordinate changes by at most b. A continuous point outside radius r thus maps only to grid paths with max |v_j^grid|>=r-b, not the original r. For nonzero q, terminal resource bins must cover the real target resource and this rounding tolerance.

Write delta_j=v_j^grid-v_j, eta_j=f_j(v_j^grid)-f_j(v_j). E delta_j=0, E delta_j^2<=b^2/4, E eta_j^2<=sigma^2 b^2/4. On a cell without a kink E eta_j=0. On a cell containing one kink, |E eta_j|<=(sigma-t_j)b/4 by the interpolation error of a hinge. For b<2B/(sigma-min t), there is at most one kink per cell. For each edge, expansion and Cauchy-Schwarz give the upper interpolation error
E E_j(v_j^grid,v_{j+1}^grid)-E_j(v_j,v_{j+1}) <= [|v_{j+1}-v_j| (|E eta_{j+1}|+|E eta_j|)+sigma b^2]/d.
No independence is assumed; the last term follows from the L2 triangle inequality on the two differences. For a site cost with piecewise C^2 bound M on its second derivative and derivative jump magnitude <=D at the kink, its expectation error is <=M b^2/8+D b/4, using linear interpolation error and subtracting a hinge to make the derivative continuous. At fixed q the analytic Gaussian site terms meet these hypotheses on the compact box, with constants yet to be validated. Add the errors to obtain a uniform E_round: if every rounded path is among the DP paths with appropriate resource band and weakened flag, the DP lower bound minus E_round is a lower bound on the continuous objective. The kinetic error alone in this crude bound at N=32,b=.025 is .265312, already greater than the .102842 relaxed-grid lower bound even before weakening the flag, so this particular mesh/bound does not certify the middle. Need sharper energy-adapted rounding, explicit q panels and tighter coordinate boxes; the exact Monge acceleration makes finer grids more feasible but does not guarantee tractability. No rounding certificate executed here.

## 4. Prior art and reproducibility

Checked abstracts of Aggarwal, Klawe, Moran, Shor, Wilber, Algorithmica 2 (1987), doi:10.1007/BF01840359; Klawe and Kleitman, SIAM J Discrete Math 3 (1990), doi:10.1137/0403009; Bein, Larmore, Park, The d-edge shortest-path problem for a Monge graph (1992), https://digital.library.unt.edu/ark:/67531/metadc1317498/. Monge searches and constrained DP are established. The specific kinetic quadrangle identity and application to a saddle selection certificate are the contributions; no novelty claim. Earlier partial averaging and Picard-Lefschetz prior-art audits still apply.

Run from experiments/2026-11-24-nova:
    OPENBLAS_NUM_THREADS=1 python structured_monge.py
    python -m py_compile structured_monge.py

Both executed successfully; NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0, Numba 0.65.1. Local copies of the imported setup and root_ball functions were used, with identical numerical formulas to the parent branch; the local setup omitted the unused audit routines and threadpoolctl call (OPENBLAS_NUM_THREADS=1 was set). Source commit 7802c382e4469e05b8c27547d4f2ea25fda2805b; verified Git blob 54c1a16f28338d6984ad153d774139b7d46cc940 matches executed script. No Lean run.

Sparks: please review the quadrangle identity, resource indexing and correlated-rounding bridge independently. In particular, try to reduce the rounding budget or find a counterexample to the proposed bounds; your independent contour work is unaffected. Next checkpoint: q-dependent energy bounds and a validated interpolation correction or exact cell-based outer relaxation. Do not claim global positivity from the lattice results. No merge.
