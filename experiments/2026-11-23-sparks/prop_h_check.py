"""Numerical test of Nova's Proposition H (real critical points, symmetric Gaussian well, a=b=0).
S(x) = x^T A x/2 + g*sum exp(-x_j^2/(2 sigma^2)), gamma = g/sigma^2.  Floating point evidence only."""
import numpy as np
from scipy.optimize import fsolve
rng = np.random.default_rng(1)
m = w = 1.0; sig = 1.0
def run(N, T, frac):
    d = T/N; n = N-1
    A = (m/d)*(np.diag((2-w*w*d*d)*np.ones(n)) - np.diag(np.ones(n-1), 1) - np.diag(np.ones(n-1), -1))
    lam = np.sort(np.linalg.eigvalsh(A)); l1 = lam[0]; l2 = lam[1] if n > 1 else np.inf
    gam = frac*l1 if frac <= 1 else (l1 + (frac-1)*((l2 if np.isfinite(l2) else 3*l1) - l1))
    g = gam*sig**2
    grad = lambda x: A@x - gam*x*np.exp(-x*x/(2*sig**2))
    hess = lambda x: A - np.diag(gam*(1-x*x/sig**2)*np.exp(-x*x/(2*sig**2)))
    roots = []
    for _ in range(3000):
        x0 = rng.normal(size=n)*rng.choice([0.3, 1.5, 4.0])
        x, info, ier, _ = fsolve(grad, x0, fprime=hess, full_output=True, xtol=1e-13)
        if ier == 1 and np.linalg.norm(grad(x)) < 1e-9 and all(np.linalg.norm(x-r) > 1e-6 for r in roots):
            roots.append(x)
    desc = []
    for r in sorted(roots, key=lambda r: r.sum()):
        ev = np.linalg.eigvalsh(hess(r)); desc.append((np.round(r, 3).tolist(), int((ev < -1e-9).sum())))
    return l1, l2, gam, desc
for (N, T) in ((2, 1.0), (3, 1.0), (4, 1.0), (6, 2.0), (9, 2.5)):
    for frac in (0.8, 1.5):   # 0.8*lambda1 (one saddle expected); halfway between lambda1 and lambda2 (three expected)
        l1, l2, gam, desc = run(N, T, frac)
        signs = all((np.array(r) > 0).all() or (np.array(r) < 0).all() or np.allclose(r, 0) for r, _ in desc)
        print(f"N={N} T={T} lam1={l1:.3f} lam2={l2:.3f} gamma={gam:.3f}: {len(desc)} real critical points, Morse indices {[i for _, i in desc]}, one-signed: {signs}")
