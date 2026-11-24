"""(a) Smooth slice-dependent tilt t(s) = a + b s + c s^2 (no separate endpoint values; a jump in t costs a term of size jump/d).
(b) For the third target at N = 8 and 16: direct multistart minimisation of h(u) = Im Phi(z_s + C u) - Im Phi(z_s) over R^(N+1),
to separate 'the path is not tilt-selected' from 'the scalar certificate is too crude'.  Floating point."""
import numpy as np
from scipy.optimize import minimize
from nslice import g, Ttot, W, W1, setup, Phi, hess, grad, reference
from nslice2 import margin, lead_oriented, saddle
def tilt(par, N):
    s = np.linspace(0, 1, N+1); return np.clip(par[0] + par[1]*s + par[2]*s*s, -0.9, 0.9)
store = {}
for label, dq, dp in (("classical target", 0, 0), ("displaced (+0.4, +0.4)", 0.4, 0.4), ("displaced (-0.6, +0.3)", -0.6, 0.3)):
    z, qb, pb = saddle(32, dq, dp); best = None
    for start in ([0.2, 0, 0], [0.1, 0.3, 0], [0.4, -0.3, 0], [0.05, 1.0, -1.0], [0.5, 0, -0.4]):
        r = minimize(lambda par: -margin(z, 32, tilt(par, 32)), start, method='Nelder-Mead', options=dict(maxiter=300, xatol=1e-3, fatol=1e-5))
        if best is None or r.fun < best.fun: best = r
    store[label] = best.x
    print(f"\n{label}: smooth tilt (a, b, c) = {np.round(best.x, 3)}")
    print("   N | margin | oriented branch valid | signed rel err hbar=.05 | hbar=.025")
    for N in (8, 16, 32, 64, 128):
        z, qb, pb = saddle(N, dq, dp); t = tilt(best.x, N); errs = []
        for hb in (0.05, 0.025):
            L, pd = lead_oriented(z, N, qb, pb, hb, t); errs.append(abs(reference(N, qb, pb, hb)/L - 1))
        print(f" {N:4d} | {margin(z, N, t):+.4f} | {pd} | {errs[0]:.4f} | {errs[1]:.4f}")
print("\nDirect test of tilt selection for the third target")
rng = np.random.default_rng(1)
for N in (8, 16):
    z, qb, pb = saddle(N, -0.6, 0.3)
    for name, t in (("uniform t = 0.2", np.full(N+1, 0.2)), ("uniform t = 0.4", np.full(N+1, 0.4)), ("optimised smooth tilt", tilt(store["displaced (-0.6, +0.3)"], N))):
        c = 1+1j*t; h0 = Phi(z, N, qb, pb).imag
        h = lambda u: Phi(z + c*u, N, qb, pb).imag - h0
        dh = lambda u: (c*grad(z + c*u, N, qb, pb)).imag
        lo = np.inf; arg = None
        for _ in range(600):
            u0 = rng.normal(0, rng.choice([0.5, 1.5, 4.0]), N+1); r = minimize(h, u0, jac=dh, method='L-BFGS-B')
            if r.fun < lo: lo, arg = r.fun, r.x
        ev = np.linalg.eigvalsh((c[:, None]*hess(z, N)*c[None, :]).imag)[0]
        print(f"   N = {N:2d}, {name}: smallest h found = {lo:+.5f} at |u| = {np.linalg.norm(arg):.3f}; smallest eigenvalue of the Hessian of h at 0 = {ev:+.4f}")
