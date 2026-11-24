"""Low-rank version of the clipped-thimble surface (Nova 037, Section 4), floating point.
z = x + i [ y_s + t o xi + sum_{k<=r} a_k clip_B(b_k . xi) ],  xi = x - x_s,  with sum a_k b_k^T the best rank-r approximation
of G = S - diag t in the mass-weighted sense (SVD of M^(1/2) G M^(-1/2)), b_k normalised so that b_k . xi is an M-coordinate.
Reports the mass-normalised local margin of the rank-r plane and the multistart minimum of h."""
import numpy as np
from scipy.linalg import sqrtm, inv, eigh
from scipy.optimize import minimize
from nslice import Phi, grad, hess, setup
from nslice2 import saddle
def build(N, r, B):
    z, qb, pb = saddle(N, -0.6, 0.3); d, w, K = setup(N); m = d*w; s = np.linspace(0, 1, N+1); t = 0.02 + 1.4*s*(1-s)
    H = hess(z, N); R = inv(sqrtm(-1j*H)); S = R.imag@inv(R.real); G = S - np.diag(t)
    U, sv, Vt = np.linalg.svd(np.sqrt(m)[:, None]*G/np.sqrt(m)[None, :])
    A = (U[:, :r]*sv[:r])/np.sqrt(m)[:, None]; Bm = Vt[:r]*np.sqrt(m)[None, :]          # G_r = A @ Bm
    C = np.eye(N+1) + 1j*(np.diag(t) + A@Bm); Mx = (C.T@H@C).imag; loc = eigh((Mx+Mx.T)/2, np.diag(m), eigvals_only=True)[0]
    xs, ys, h0 = z.real, z.imag, Phi(z, N, qb, pb).imag
    def zmap(x):
        q = Bm@(x-xs); return x + 1j*(ys + t*(x-xs) + A@np.clip(q, -B, B)), (abs(q) < B)
    def h(x): return Phi(zmap(x)[0], N, qb, pb).imag - h0
    def dh(x):
        zz, mask = zmap(x); J = np.diag(t) + A@(mask[:, None]*Bm); gr = grad(zz, N, qb, pb); return gr.imag + J.T@gr.real
    return h, dh, xs, loc, sv[:6]
for N in (32, 64):
    print(f"N = {N}: leading singular values of the mass-weighted S - diag t: {np.round(build(N, 1, 0.1)[4], 3)}")
    for r in (1, 2, 3, 5):
        for B in (0.1, 0.3, 1.0):
            h, dh, xs, loc, sv = build(N, r, B); rng = np.random.default_rng(11); lo = np.inf; arg = xs
            for k in range(300 if N == 32 else 150):
                x0 = xs + rng.normal(0, rng.choice([0.2, 0.7, 2.0, 5.0]), N+1); res = minimize(h, x0, jac=dh, method='L-BFGS-B')
                if res.fun < lo: lo, arg = res.fun, res.x
            print(f"   rank {r}, B = {B}: local margin {loc:+.4f};  smallest h found {lo:+.6f} at distance {np.linalg.norm(arg-xs):.3f}")
