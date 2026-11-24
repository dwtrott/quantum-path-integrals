"""Smooth low-mode version of the adapted surface, floating point.  Replace G = S - diag t by its compression to the r
smoothest modes: G_r = P G P with P the M-orthogonal projector on the first r generalized eigenvectors of the background
imaginary quadratic Hessian B0 = (1 - t^2) D + K t + t K (no kicks) relative to M.  Surface:
z = x + i [ y_s + t o xi + E clip_B( F xi ) ],  with G_r = E F,  F xi = mode coordinates of G_r xi  (r numbers).
Failed predecessor: low_rank.py (SVD truncation of G), whose rank-r planes are rough in time and have local margins below -100."""
import numpy as np
from scipy.linalg import sqrtm, inv, eigh
from scipy.optimize import minimize
from nslice import Phi, grad, hess, setup
from nslice2 import saddle
def build(N, r, B):
    z, qb, pb = saddle(N, -0.6, 0.3); d, w, K = setup(N); m = d*w; s = np.linspace(0, 1, N+1); t = 0.02 + 1.4*s*(1-s)
    H = hess(z, N); R = inv(sqrtm(-1j*H)); S = R.imag@inv(R.real); G = S - np.diag(t)
    Dm = np.zeros((N+1, N+1)); Dm[0, 0] = Dm[-1, -1] = 1; T = np.diag(t); B0 = Dm + K@T + T@K - T@Dm@T
    lam, V = eigh((B0+B0.T)/2, np.diag(m)); V = V[:, :r]                      # V^T M V = I
    E = V; F = (V.T*m[None, :])@G@V@(V.T*m[None, :])                          # G_r = V (V^T M G V) V^T M
    C = np.eye(N+1) + 1j*(T + E@F); Mx = (C.T@H@C).imag; loc = eigh((Mx+Mx.T)/2, np.diag(m), eigvals_only=True)[0]
    xs, ys, h0 = z.real, z.imag, Phi(z, N, qb, pb).imag
    def zmap(x):
        q = F@(x-xs); return x + 1j*(ys + t*(x-xs) + E@np.clip(q, -B, B)), (abs(q) < B)
    def h(x): return Phi(zmap(x)[0], N, qb, pb).imag - h0
    def dh(x):
        zz, mask = zmap(x); J = T + E@(mask[:, None]*F); gr = grad(zz, N, qb, pb); return gr.imag + J.T@gr.real
    return h, dh, xs, loc
for N in (32, 64):
    print(f"N = {N}")
    for r in (1, 2, 3, 5, 8):
        loc = build(N, r, 0.1)[3]; line = f"   r = {r}: local margin {loc:+.4f};"
        if loc > 0:
            for B in (0.1, 0.3, 1.0):
                h, dh, xs, _ = build(N, r, B); rng = np.random.default_rng(11); lo = np.inf; arg = xs
                for k in range(300 if N == 32 else 150):
                    x0 = xs + rng.normal(0, rng.choice([0.2, 0.7, 2.0, 5.0]), N+1); res = minimize(h, x0, jac=dh, method='L-BFGS-B')
                    if res.fun < lo: lo, arg = res.fun, res.x
                line += f"  B = {B}: min h {lo:+.5f} (dist {np.linalg.norm(arg-xs):.2f});"
        print(line)
