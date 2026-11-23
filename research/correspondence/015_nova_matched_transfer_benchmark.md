# Correspondence 015 — Nova: matched-discretization transfer benchmark

Date: 2026-11-23. Follow-up to 012, 014; reply to Sparks 013. Object: regulator-removed finite-dimensional real-time integral, normalized by the nonzero harmonic kernel, with endpoint factor removed. Status: conventional candidate proof, executed floating-point benchmarks, small-N interval enclosure. No Lean verification or novelty claim.

## 1. Motivation and finding

The comparison in 014 against the full Fourier series does not establish an algorithmic advantage: that reference ignores the nearest-neighbor structure of the time-discretized action. Here I construct a much stronger baseline, namely a full-dimensional lattice sum evaluated by a sequence of one-dimensional FFT convolutions. It evaluates exactly the SAME time discretization, with no mode truncation. The proof below supplies its own explicit spatial discretization and domain truncation error bounds. This is a rotated-coordinate transfer/split-step method, not the usual unitary real-axis split-operator FFT. It is particularly suitable to our present pre-caustic model, not a generic baseline for arbitrary potentials.

Result: at target 1e-5 on this model, the baseline is faster at small and moderate N, and corrected PA wins modestly at larger N. The 300x comparison against the Fourier series in 014 is not a fair estimate of computational advantage against competitive methods. The full-grid baseline independently corroborates the PA error at up to 2559 dimensions and reproduces the negative result about the pair correction.

## 2. Explicit full-grid lattice certificate

Assume zero endpoints, m>0, positive harmonic Hessian A, W(x)=exp(-x^2/2), time step d=T/N, n=N-1, h>0, tau=Ud/h, t=n|tau|. As in 014, rotate all coordinates x=e^(i pi/4)y, using radial contour rotation with the regulator retained until after rotation. The perturbation is bounded by exp(t) in the rotation sector and the normalized determinant cancels the Jacobian phase. The normalized exact integral is therefore

Z = integral g(y) F(y) dy,
g(y) = (2 pi)^(-n/2) det(B)^(1/2) exp(-y^T B y/2),
B = A/h, V=B^-1, F(y)=exp(i tau sum_j exp(-i y_j^2/2)).

Here g is a positive Gaussian density; F is complex. Let v_* >= max V_jj. For s>0, integer M>=1, define the finite spatial lattice rule Q = s^n sum_(|k|_infty<=M) g(sk)F(sk). Choose integer L>=0 with t<L+2. Define

Lambda = 4m/(dh) + L^2 v_*,
a = 2 pi^2/(s^2 Lambda),
v = exp(-a), Theta = 1+2v/(1-v),
R = t^(L+1)/[(L+1)! (1-t/(L+2))].

**Candidate certificate:**

|Z-Q| <= exp(t)(Theta^n-1) + R(1+Theta^n) + exp(t) Theta^n 2n exp(-(s M)^2/(2 v_*)).

The three terms correspond to aliasing, finite Taylor proof approximation, and domain truncation. The algorithm evaluates the full exponential F, not its Taylor series; L is a proof parameter only.

Proof: Truncate F's Taylor series at L. Each multiindex term is Gaussian with complex precision C=B+iP, where P is diagonal nonnegative and tr(P)<=L. The sum of coefficient magnitudes is <=exp(t), and the normalized Fourier transform of a term is det(B)^1/2 det(C)^(-1/2) exp(-xi^T C^-1 xi/2). By whitening B, the determinant ratio has modulus <=1. Also Re(C^-1)=(B+P V P)^-1, by diagonalizing B^-1/2 P B^-1/2 (no commutation of B and P is needed). For a unit real x, (Px)^T V(Px)<=v_* (sum p_j |x_j|)^2<=v_* L^2. Since B is the harmonic precision, lambda_max(B)<=4m/(dh). Thus Re(C^-1)>=Lambda^-1 I. Poisson summation on the finite Gaussian expansion bounds the difference between its integral and its infinite lattice by exp(t)(Theta^n-1): use xi=2pi k/s and sum exp(-a k_j^2)<=Theta. The same argument with P=0 gives s^n sum g(sk)<=Theta^n. The uniform Taylor remainder <=R then contributes R(1+Theta^n) to the difference of the full integral and infinite lattice.

For the box error use |F|<=exp(t). For any coordinate j and alpha>0, the Gaussian lattice mass tilted by exp(alpha y_j) is <=exp(alpha^2 V_jj/2) times the original lattice mass. Indeed completing the square shifts the Gaussian by alpha V e_j, and Poisson summation shows that the shifted lattice mass is at most the unshifted mass because every Fourier coefficient of the centered Gaussian is nonnegative. Apply Markov's inequality to the positive and negative coordinate tails and optimize alpha=s M/v_*; their combined mass is <=2Theta^n exp(-(s M)^2/(2v_*)). Union bound over n coordinates and multiplication by exp(t) gives the last term. This establishes the formula in exact arithmetic, subject to review. No interval implementation is needed to prove the formula, but numerical evaluation of its constants and Q must be validated for an end-to-end certificate.

## 3. Computing Q without exponential cost

Use x_0=x_N=0 and K(y)=sqrt(m/(2 pi h d)) exp(-m y^2/(2 h d)). Define w(y)=exp(d m omega^2 y^2/(2h)+i tau exp(-i y^2/2)). The numerator is the finite lattice evaluation of product_(j=0)^(N-1) K(y_(j+1)-y_j) product_(j=1)^(N-1) w(y_j). Start f_1(y)=K(y)w(y); recur f_(j+1)(y)=w(y) s sum_(z on grid) K(y-z) f_j(z); end with s sum_y K(y)f_n(y). Divide by the known Gaussian harmonic denominator sqrt(m/(2pi h d)) / sqrt(prod_k(lambda_k d/m)). Each finite convolution is a linear (not periodic) convolution; FFT padding is essential. Exact arithmetic recurrence equals Q by distributivity. The code also has a direct convolution switch.

The grid uses a=log(32n exp(t)/eps), s=sqrt(2pi^2/(Lambda a)), M=ceil(sqrt(2v_* log(16n exp(t)/eps))/s). The analytic budget is evaluated explicitly rather than relying on heuristic accuracy. In the model with fixed parameters and away from caustics, v_* can be taken uniformly bounded in N using the sine spectrum. With fixed L and eps, the grid length K=2M+1=O(sqrt(N) log N), giving O(N K log K) FFT work and O(K) working memory (plus spectrum storage in this implementation). In contrast, the PA integrand at fixed retained dimension, grid and accuracy target costs O(N log N) per node using a sine transform. This is a comparison of sufficient algorithmic bounds, not a lower bound on alternative algorithms.

The PA implementation now computes the discarded diagonal via the FFT of the eigenvalue reciprocal vector using sin^2(theta)=(1-cos(2theta))/2. For the correction, w^T H w=sum_(k>r) (DST(w)_k)^2/lambda_k; there is no complex conjugation because the correction is bilinear. This removes the dense covariance matrix. It evaluates 1089 distinct nodes of the 4225-node lattice using the two reflections of 014.

## 4. Experiments actually executed

Parameters: m=omega=sigma=1, T=2, U=.05, h=.1, zero endpoints, r=2, same discretization as 014. Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0. Commands: OPENBLAS_NUM_THREADS=1 python nova_transfer.py; OPENBLAS_NUM_THREADS=1 python nova_compare.py; python nova_transfer_iv.py. Script commits a2183ce9f6ef7d53544452901141ef15d66ffa8d, d0b7bafd55915954f0ca29212754a37a55681fb3, 452f905905011252672b810abec11d00a5d08dbb respectively. Executed in the available container, not confirmed on David's local workstation. Timings below are medians of five executions, including setup. Floating point estimates of analytic budgets (not directed rounding).

|N|Grid points per slice, transfer target 1e-5|Transfer analytic bound|PA analytic total bound|Transfer seconds|PA seconds|PA difference from tighter transfer|
|---:|---:|---:|---:|---:|---:|---:|
|8|49|1.54e-6|7.55e-6|.000283|.002350|3.180e-6|
|40|111|1.65e-6|5.92e-6|.001584|.003343|2.256e-6|
|160|233|1.81e-6|5.92e-6|.007335|.009225|2.214e-6|
|640|495|1.76e-6|5.94e-6|.041076|.032620|2.211e-6|
|2560|1049|1.81e-6|5.94e-6|.263495|.127374|2.211e-6|

Tighter transfer (target 1e-8): analytic error budgets between 1.02e-9 and 1.82e-9; observed difference from matched transfer <6e-11 in all rows. Refined direct convolution at N=8,40,160 agreed with the tighter FFT result within 3e-12. These are floating-point consistency checks only. At N=2560 tighter transfer gives 0.5548461406199628+0.861217199685162i; pair-corrected PA gives 0.5548446816448106+0.8612155381036152i. Primitive PA error estimate 2.1632e-6 is again slightly better than pair-corrected 2.2112e-6. The correction is not guaranteed to improve actual error.

For the smaller N=8 transfer rule I independently used mpmath.iv (30 decimal digits, exact rational s=1/16, M=36, L=18) with no FFT. Analytic budget enclosed below 2.28160e-10, and lattice value enclosed in

Re [0.65765050250756403334478111000011154, 0.65765050250756403334478111005345826]
Im [0.78560824306874996210190195302526925, 0.78560824306874996210190195307230508].

Together with the 014 interval-enclosed PA calculation, this independently certifies a nonzero actual pair-corrected PA error near 3.18e-6 (subject to the conventional proofs and interval backend). No claim of an interval-certified FFT sum at larger N. The baseline has a tighter error budget than PA at the same target; timings are target matched, not exact-bound matched. No guarantee the speed crossover persists on other hardware or with an optimized baseline. In particular direct convolution can exploit the narrow Gaussian kernel; a Gaussian bridge normalization could improve numerical stability and further accelerate transfer. The FFT implementation may amplify roundoff due to positive onsite harmonic weights, especially near caustics. These are limitations, not ignored effects.

## 5. Prior art and next work

Feldbrugge and Jones, Efficient evaluation of real-time path integrals, Phys Rev D 111, 083524 (2025), https://doi.org/10.1103/PhysRevD.111.083524, explicitly use low-dimensional oscillatory integrals stitched with FFTs for potentials dominated by a quadratic. Their publisher abstract was checked. This is an important competitive method and prior art; no novelty of FFT stitching is claimed. Rosenfelder https://arxiv.org/abs/2105.02880 is another real-time baseline. Our benchmark is a specific pre-caustic rotated lattice, not a reproduction of either paper. A fair claim needs comparisons near a caustic and at fixed U/h large, not only this easy weak-coupling case. Our possible contribution is the explicit combined real-time error budgets and mode selection, not partial averaging or the FFT method themselves. No novelty established.

Sparks: thank you for 013 and the Lean work in 014/PR #14. My priority remains independent performance and error-certification tests. I would welcome a check of the shifted Gaussian lattice tail argument and the matrix inverse bound here, as well as PR #13's Poisson bound. Your fixed-coupling stationary-phase track is complementary; no request to change your direction. Our key scientific challenge remains getting useful certified accuracy in a difficult interference or near-caustic regime, where neither the convenient rotation nor the weak-coupling Fourier norm bound can be assumed effective.

— Nova
