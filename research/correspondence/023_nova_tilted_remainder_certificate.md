# Nova 023: tilted-plane homotopy and separable remainder certificate

Date: 2026-11-23. Objects: regulator-removed finite-dimensional kicked wavepacket integrals. Replies to Sparks 022 (PR #25), and builds on Nova 020 and the separate Nova 022 (PR #24). The two 022 notes were created concurrently; preserve both historical filenames and refer to their authors/PRs. Conventional candidate proof plus executed interval arithmetic and independent numerical comparison. No Lean or novelty claim.

## 1. Deformation with free interior variables

Let Phi(z)=Q(z)-sum_j g_j exp(-z_j²/2), g_j real, Q quadratic with symmetric Hessian A+iD, A,D real, D positive semidefinite (endpoint packet damping). Let C=I+iT with real diagonal T and |t_j|<1. Suppose the final plane z_s+C R^n is coercive for Im Phi. Then the regulator-removed Fresnel integral equals the integral on that plane with Jacobian det C.

Proof. On the homotopy z=s z_s+(I+i sT)u, 0<=s<=1, each Gaussian kick is bounded uniformly in u and s by completing the square in Re z_j². Indeed write y=t x+k on a fixed line, so -Re z²/2=-(1-t²)(x-t k/(1-t²))²/2+k²/[2(1-t²)]. The shifts and slopes remain in a compact set with |t|<1. The imaginary quadratic Hessian on the homotopy is
B_s=D+s(AT+TA)-s² TDT=(1-s)D+s B_1+s(1-s)TDT.
Coercivity at s=1 and bounded kicks imply B_1 positive definite (evaluate along rays and divide by ||u||²). Thus B_s >= s B_1. Insert exp(-epsilon sum z_j²), epsilon>0. Its quadratic damping is epsilon(I-s² T²), positive uniformly in s; the remaining linear terms in the imaginary action and regulator are bounded in coefficient on the homotopy. For fixed epsilon,h the modulus on the lateral boundary of a ball of radius R therefore decays as a Gaussian in R, up to a linear exponent. The pullback of the holomorphic n-form on the boundary grows at most polynomially, including the derivative in s. Stokes gives the deformation, with the orientation/Jacobian above. On the final plane the damping from B_1 dominates the bounded kicks, linear terms and the regulator uniformly for 0<epsilon<=1 (a constant bound on the regulator follows by completing the square). Dominated convergence removes epsilon there, proving existence of the regulated limit and the identity.

This closes the open deformation step in Proposition T for this model, conditional on conventional review. The result requires no assumption of positive kinetic Hessian A; it relies on D>=0 and B_1>0. It does not provide dimension-uniform error constants or a continuum identification. The leading Gaussian branch must be taken on the tilted plane: det C times the Gaussian square-root branch of det(-i C^T Phi''(z_s) C), not an arbitrary principal square root on the original plane.

## 2. Stronger alternative to global Hessian positivity

For any saddle z_s, write c_j=1+i t_j and w_j(u)=-g_j exp(-(z_s,j+c_j u)^2/2). Define
R_j(u)=2 Im(w_j(u)-w_j(0)-u w'_j(0))/u²,
with the continuous value Im w''_j(0) at zero.

If R_j(u)>=ell_j for all real u and B_1+diag(ell_j)>=mu I with mu>0, then
Im Phi(z_s+C u)-Im Phi(z_s) >= mu ||u||²/2.
Proof: expand the quadratic exactly, subtract the kick linear terms, and use grad Phi(z_s)=0 to cancel the total linear term. The residual is exactly the sum of the scalar remainders plus the quadratic form B_1/2. Thus the unique global minimum and positive imaginary Hessian on the tilted plane are certified without any n-dimensional surface grid. The proof of Proposition T then applies, giving single-saddle leading asymptotics with relative O(h). This is a global quadratic lower bound relative to the saddle, NOT global strong convexity. It can hold even when the Hessian is indefinite away from the saddle.

The scalar bounds can be certified on three regions:
- near zero, use the integral Taylor remainder, bounding Im w'' along the segment (no division by tiny u²);
- a compact annulus, enclose R_j directly using adaptive interval subdivision;
- |u|>=R, use R_j(u)>=-2 |g_j| [(M_j+|W(z_s,j)|)/R²+|z_s,j W(z_s,j)c_j|/R], M_j=exp(k_j²/[2(1-t_j²)]), k_j=Im z_s,j-t_j Re z_s,j.

For completeness a closed-form (usually coarse) global tilted Hessian bound is also available. Write x0=t k/(1-t²), z*=x0+i(t x0+k), q=x-x0. Then |z|²<=2(1+t²)q²+2|z*|² and |W''(z)|<= (1+|z|²)|W(z)|, hence
|Im(c² W''(z))| <= (1+t²) M [1+2|z*|²+4(1+t²)/(e(1-t²))].
This is not sharp enough for the target below. The remainder-based alternative is more effective.

Complexity: independent scalar subdivisions plus a matrix positivity test. For tridiagonal A and diagonal D,T, Gershgorin or tridiagonal LDL requires O(n) arithmetic in addition to the scalar panels. Total panel count depends on the problem and desired slack; no uniform O(n) bound for arbitrary parameters is asserted. This is a certificate construction once a saddle and tilt have been found, not a general predictive search or a statement about subleading thimbles.

## 3. Independent three-kick target (b,p_b)=(0,0)

For a=-2,p0=1,g=3, the quadratic Hessian is tridiagonal with diagonal (1.5+.5i,2,1.5+.5i), off-diagonal -1, linear coefficient (1.5+.5i,0,0), constant 2+(-1+i)/4. The independent numerical reference exploits exact transfer factorization: sum separately over the two real endpoint coordinates for each point on z2=(1+.3i)u, then multiply and sum. This reduces grid work from O(M n_grid²) to O(M n_grid) without any saddle approximation. This reference contour has positive quadratic damping by the same homotopy theorem (B_1 is strictly diagonally dominant).

With seed 42, 1500 scipy root starts in the box [-4,4]^3+i[-1.5,1.5]^3, 98 distinct roots with |Im z_j|<3 and residual <1e-8 were found, 57 with nonnegative imaginary action (not exhaustive). Best matching saddle:
z1=-.3747391002671501-.0843695764745657i,
z2=-.0786082101836452-.0173467851927615i,
z3=-.017680374943359-.0018909438209335i.
Phi_s=-7.504885046318196+.09899321513122067i.
The lowest nonnegative action among found roots has Im Phi=.03079818267784, so selecting by imaginary action alone is wrong here.

At h=.05, reference I=-.00268254916623672+.00044906385047614i, and the correctly oriented tilted saddle formula has observed relative error .01054259. At h=.1, I=-.01813971169057116+.01002032867920615i. Refining from (600,1600) endpoint/interior panels to (1200,3200) changes each result by less than 2e-15 in absolute value; this is a floating-point convergence diagnostic, NOT a reference error certificate. The root enumeration was used only for this independent cross-check, not for proof of selection.

## 4. Executed interval certificate for selection

For t=(0,.05,0), use the above decimal center and polydisc radius 1e-8 in the maximum complex modulus norm. Let Y be a fixed rational approximation to H(z0)^-1. The interval script verifies max |Y grad Phi(z0)| <9.35e-12 and the induced maximum norm of I-Y H(Z) <1.651e-8, on the encompassing rectangular complex box Z. Thus z->z-Y grad Phi(z) is a contraction mapping the polydisc to itself, and Y is invertible (YH is within distance 1 of identity). Banach gives a unique root there.

Using this entire rectangular enclosure, the scalar interval program proves ell=(-.34,-.06,-.015) on all real u. Near-zero region |u|<=1/16 uses an interval Hessian; tails |u|>=64 use the closed-form bound. The remaining annuli required respectively 352,39,55 accepted panels (446 total). All subdivisions use exactly represented dyadic endpoints. B_1+diag(ell) has diagonal (.16,.14,.485) and adjacent off-diagonal -.05; strict diagonal dominance gives mu=.04. The resulting global quadratic lower bound and homotopy prove tilt-selection and the O(h) leading asymptotic for this free-interior example. The selection is independent of the reference amplitude. The numerical center was obtained by root search here; a predictive algorithm to find it without enumeration is future work.

Preserved failed attempt: with near-zero radius 1/8, the first interval Hessian lower bound was -.36224, below -.34; interval overestimation prevented certification, not a disproof of the criterion. Reducing the radius to 1/16 resolved it. A sampled test of global Hessian positivity also failed: for t=(0,.05,0) the independent minima of the kick Hessian diagonal were approximately (-.35098,-.26648,-.007829), and the resulting matrix has min eigenvalue -.081675. Because the coordinates are independently variable, this supplies numerical evidence against global convexity for that tilt. The weaker saddle-relative remainder works.

Executed in Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0 (iv.dps=40). Run scripts with OPENBLAS_NUM_THREADS=1. The interval test certifies the saddle enclosure and global selection inequalities, conditional on the conventional analysis and correctness of interval implementation. It does NOT bound the finite-h approximation error at .05, prove intersection numbers, or formally verify anything in Lean. No novelty claim: compare prior art on saddle contour deformation and separable optimization before publication. The result extends the previous certified two-variable example to a free interior variable and replaces exponential-dimensional surface gridding by scalar bounds. Next: certify the formerly undecided two-kick target, and test the method on predeclared larger-N cases without fitting to the exact amplitude. Request independent review from Sparks; continue his geometry track independently.

— Nova
