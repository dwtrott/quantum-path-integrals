"""Local plane C = I + i (sigma I + V A V^T M): uniform slope sigma on every slice plus a correction on the r smoothest modes.
How much of the thimble-plane local margin (0.72) does this structured family recover?  Floating point, third target."""
import numpy as np
from scipy.linalg import eigh, sqrtm, inv
from scipy.optimize import minimize
from nslice import hess, setup, classical, newton
def saddle(N, dq, dp):
    x, p = classical(N); z = x.astype(complex)
    for s in np.linspace(0, 1, 21)[1:]: z = newton(z, N, x[-1]+s*dq, p+s*dp, it=12)
    return newton(z, N, x[-1]+dq, p+dp)
for N in (32, 64):
    z = saddle(N, -0.6, 0.3); d, w, K = setup(N); m = d*w; H = hess(z, N); I = np.eye(N+1)
    R = inv(sqrtm(-1j*H)); S = R.imag@inv(R.real); sv = np.linalg.svd(np.sqrt(m)[:, None]*S/np.sqrt(m)[None, :], compute_uv=False)
    print(f"N = {N}: mass-weighted singular values of S: largest {np.round(sv[:3], 3)}, median {np.median(sv):.3f}, smallest {np.round(sv[-3:], 3)}; diagonal of S: min {S.diagonal().min():.3f}, max {S.diagonal().max():.3f}")
    def marg(C): Mx = (C.T@H@C).imag; return eigh((Mx+Mx.T)/2, np.diag(m), eigvals_only=True)[0]
    V_all = eigh(K + np.diag(m), np.diag(m))[1]
    for r in (0, 1, 2, 3, 5):
        V = V_all[:, :r]; P = V.T*m[None, :]; best = None
        for sig0 in (0.3, 0.6, 0.9):
            f = lambda a: -marg(I + 1j*(np.clip(a[0], 0.02, 0.99)*I + (V@a[1:].reshape(r, r)@P if r else 0)))
            res = minimize(f, np.r_[sig0, np.zeros(r*r)], method='Nelder-Mead', options=dict(maxiter=6000, xatol=1e-4, fatol=1e-7))
            if best is None or res.fun < best.fun: best = res
        print(f"   r = {r}: best local margin {-best.fun:+.4f} at sigma = {np.clip(best.x[0], 0.02, 0.99):.3f}" + (f", largest |A| entry {abs(best.x[1:]).max():.2f}" if r else ""))
