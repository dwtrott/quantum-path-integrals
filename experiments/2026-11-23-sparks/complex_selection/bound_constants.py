"""Constants of the proposed N-uniform bound |I/L - 1| <= (C1 + C2) hbar + (exponentially small), classical target, floating point.
Plane z = z_s + R u, R = (-iH)^(-1/2), u ~ N(0, hbar I).  delta_j = (R u)_j, E|delta_j|^2 = hbar s_j^2, s_j^2 = (R R^*)_jj.
C1 = g M4 sum_j m_j s_j^4 / 8                       (mean of the quartic remainder; the cubic term has zero mean)
C2 = F_c g^2 (15/72) ( sum_j m_j M3 s'_j^3 )^2       (second-order remainder under the dominating Gaussian with covariance
                                                     hbar (I - c X^T M X)^(-1), X = Re R, F_c = det(I - c X^T M X)^(-1/2))
M3, M4: suprema of |W'''|, |W''''| on the strip |Im z| <= Y.  The hypothesis h(u) >= (|u|^2 - c ||X u||_M^2)/2 is NOT verified here;
the script estimates the smallest admissible c on the plane by multistart maximisation of (|u|^2 - 2 h(u)) / ||X u||_M^2."""
import numpy as np
from scipy.linalg import sqrtm, inv, eigh
from scipy.optimize import minimize
from nslice import g, Phi, grad, hess, setup, W
from order_hbar import saddle, a1, W3, W4
def sup_strip(f, Y):
    x = np.linspace(-6, 6, 2401); y = np.linspace(-Y, Y, 41); X, Yy = np.meshgrid(x, y); return np.abs(f(X + 1j*Yy)).max()
for Y in (0.0, 0.2):
    print(f"strip half-width Y = {Y}: M3 = {sup_strip(W3, Y):.4f}, M4 = {sup_strip(W4, Y):.4f}")
M3, M4 = sup_strip(W3, 0.2), sup_strip(W4, 0.2)
print("    N | max s_j^2 | C1 | F_c, C2 at c = 0.5 | F_c, C2 at c = 1.0 | |a_1|")
for N in (16, 32, 64, 128, 256):
    z, qb, pb = saddle(N, 0, 0); H = hess(z, N); d, w, K = setup(N); m = d*w; R = inv(sqrtm(-1j*H)); X = R.real
    s2 = np.sum(abs(R)**2, axis=1); C1 = g*M4*np.sum(m*s2**2)/8; out = []
    for c in (0.5, 1.0):
        A = np.eye(N+1) - c*X.T@np.diag(m)@X; Fc = np.linalg.det(A)**-0.5; Ai = inv(A); sp2 = np.einsum('jk,kl,jl->j', R, Ai, R.conj()).real
        out.append(f"{Fc:.4f}, {Fc*g*g*15/72*np.sum(m*M3*sp2**1.5)**2:.3f}")
    print(f"  {N:4d} | {s2.max():.4f} | {C1:.4f} | " + " | ".join(out) + f" | {abs(a1(z, N)):.4f}")
print("estimate of the smallest admissible c on the plane (multistart; NOT a certificate, and the plane is not the admissible surface)")
rng = np.random.default_rng(3)
for N in (16, 32):
    z, qb, pb = saddle(N, 0, 0); H = hess(z, N); d, w, K = setup(N); m = d*w; R = inv(sqrtm(-1j*H)); X = R.real; h0 = Phi(z, N, qb, pb).imag
    def negratio(u):
        hv = Phi(z + R@u, N, qb, pb).imag - h0; den = (X@u)@(m*(X@u)); return -(u@u - 2*hv)/den
    best = 0; arg = None
    for k in range(300):
        u0 = rng.normal(0, rng.choice([0.05, 0.2, 0.5, 1.0]), N+1); r = minimize(negratio, u0, method='L-BFGS-B', bounds=[(-1.5, 1.5)]*(N+1))
        if -r.fun > best: best, arg = -r.fun, r.x
    lam1 = 1/eigh(X.T@np.diag(m)@X, eigvals_only=True)[-1]
    print(f"   N = {N}: largest ratio found {best:.4f} within |u_k| <= 1.5 (at max |Im(R u)| = {np.max(abs((R@arg).imag)):.3f});  lambda_1 = {lam1:.4f}")
