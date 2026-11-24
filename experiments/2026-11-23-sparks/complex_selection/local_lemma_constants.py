"""Constants of Lemma E-local (correspondence 049), classical target, floating point.
P = X^T M X + Y^T M Y with R = X + iY = (-iH)^(-1/2);  kappa = g M3(r) r / 3;  requires kappa ||P|| < 1.
C1 = g M4(r) sum m_j s_j^4 / 8,  s_j^2 = (R R^*)_jj
C2 = det(I - kappa P)^(-1/2) (g M3(r))^2 (15/72) (sum m_j s'_j^3)^2,  s'_j^2 = (R (I - kappa P)^(-1) R^*)_jj
Tail by the union bound: 1 - P(S_r) <= 4 (N+1) exp(-r^2 / (4 hbar s_max^2))."""
import numpy as np
from scipy.linalg import sqrtm, inv
from nslice import g, hess, setup
from order_hbar import saddle, a1, W3, W4
def sup_strip(f, Y):
    x = np.linspace(-6, 6, 2401); y = np.linspace(-Y, Y, 41); X, Yy = np.meshgrid(x, y); return np.abs(f(X + 1j*Yy)).max()
print("    N |  r   | ||P|| | tr P | kappa ||P|| | Fredholm factor | C1 | C2 | C1 + C2 | |a_1| | tail at hbar = 1/40, 1/160")
for N in (32, 64, 128, 256):
    z, qb, pb = saddle(N, 0, 0); H = hess(z, N); d, w, K = setup(N); m = d*w; R = inv(sqrtm(-1j*H)); X, Y = R.real, R.imag
    P = X.T@np.diag(m)@X + Y.T@np.diag(m)@Y; nP = np.linalg.eigvalsh(P)[-1]; s2 = np.sum(abs(R)**2, axis=1)
    for r in (0.1, 0.2, 0.3):
        M3, M4 = sup_strip(W3, r), sup_strip(W4, r); kap = g*M3*r/3
        if kap*nP >= 1: print(f"  {N:4d} | {r} | kappa ||P|| = {kap*nP:.3f} >= 1, lemma not applicable"); continue
        A = np.eye(N+1) - kap*P; F = np.linalg.det(A)**-0.5; sp2 = np.einsum('jk,kl,jl->j', R, inv(A), R.conj()).real
        C1 = g*M4*np.sum(m*s2**2)/8; C2 = F*(g*M3)**2*15/72*np.sum(m*sp2**1.5)**2
        tails = [4*(N+1)*np.exp(-r*r/(4*hb*s2.max())) for hb in (1/40, 1/160)]
        print(f"  {N:4d} | {r} | {nP:.4f} | {np.trace(P):.4f} | {kap*nP:.3f} | {F:.4f} | {C1:.3f} | {C2:.2f} | {C1+C2:.2f} | {abs(a1(z, N)):.4f} | {tails[0]:.2e}, {tails[1]:.2e}")
