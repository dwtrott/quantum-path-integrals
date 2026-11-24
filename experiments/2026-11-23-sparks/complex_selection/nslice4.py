"""(a) Smooth tilt constrained to 0.02 <= t(s) <= 0.9, optimised at N = 64, evaluated to N = 256.
(b) Local obstruction test for the third target: maximise over all diagonal tilts (one free t_j per slice, |t_j| <= 0.95)
the smallest eigenvalue of Im(C H C), the Hessian at the saddle of the imaginary action on the tilted plane.  Floating point."""
import numpy as np
from scipy.optimize import minimize
from nslice import hess, reference
from nslice2 import margin, lead_oriented, saddle
def tilt(par, N):
    s = np.linspace(0, 1, N+1); return np.clip(par[0] + par[1]*s + par[2]*s*s, 0.02, 0.9)
for label, dq, dp in (("classical target", 0, 0), ("displaced (+0.4, +0.4)", 0.4, 0.4)):
    z, qb, pb = saddle(64, dq, dp); best = None
    for start in ([0.2, 0, 0], [0.1, 0.3, 0], [0.4, -0.3, 0], [0.05, 1.0, -1.0], [0.3, 0.4, -0.4]):
        r = minimize(lambda par: -margin(z, 64, tilt(par, 64)), start, method='Nelder-Mead', options=dict(maxiter=300, xatol=1e-3, fatol=1e-5))
        if best is None or r.fun < best.fun: best = r
    print(f"\n{label}: constrained smooth tilt (a, b, c) = {np.round(best.x, 3)}")
    print("   N | margin | oriented branch valid | signed rel err hbar=.05 | hbar=.025")
    for N in (8, 16, 32, 64, 128, 256):
        z, qb, pb = saddle(N, dq, dp); t = tilt(best.x, N); errs = []
        for hb in (0.05, 0.025):
            L, pd = lead_oriented(z, N, qb, pb, hb, t); errs.append(abs(reference(N, qb, pb, hb)/L - 1))
        print(f" {N:4d} | {margin(z, N, t):+.4f} | {pd} | {errs[0]:.4f} | {errs[1]:.4f}")
print("\nLocal obstruction test, third target (-0.6, +0.3)")
rng = np.random.default_rng(2)
for N in (8, 16, 32):
    z, qb, pb = saddle(N, -0.6, 0.3); H = hess(z, N)
    def f(t):
        c = 1+1j*np.clip(t, -0.95, 0.95); return -np.linalg.eigvalsh((c[:, None]*H*c[None, :]).imag)[0]
    best = np.inf
    for k in range(12):
        t0 = rng.uniform(-0.3, 0.6, N+1) if k else np.full(N+1, 0.2)
        r = minimize(f, t0, method='Powell', options=dict(maxiter=20000, xtol=1e-4, ftol=1e-7)); best = min(best, r.fun)
    HR = np.linalg.eigvalsh(H.real); print(f"   N = {N:2d}: largest achievable smallest eigenvalue of Im(C H C) over diagonal tilts = {-best:+.5f};  eigenvalues of Re H below zero: {int(np.sum(HR < 0))};  smallest eigenvalue of Im H = {np.linalg.eigvalsh(H.imag)[0]:+.4f}")
for label, dq, dp in (("classical target", 0, 0), ("displaced (+0.4, +0.4)", 0.4, 0.4)):
    z, qb, pb = saddle(32, dq, dp); print(f"   for comparison, {label}, N = 32: eigenvalues of Re H below zero: {int(np.sum(np.linalg.eigvalsh(hess(z, 32).real) < 0))}")
