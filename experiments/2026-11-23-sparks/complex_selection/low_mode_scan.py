"""Scaled low-mode correction: C = I + i (T + eps * V A V^T M), V the r smoothest M-orthonormal modes, A an r x r real matrix.
Maximise the mass-normalised local margin over A (free), starting from A = 0.  Floating point."""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from nslice import hess, setup
from nslice2 import saddle
for N in (32, 64, 128):
    z, qb, pb = saddle(N, -0.6, 0.3); d, w, K = setup(N); m = d*w; s = np.linspace(0, 1, N+1); t = 0.02 + 1.4*s*(1-s); H = hess(z, N)
    Dm = np.zeros((N+1, N+1)); Dm[0, 0] = Dm[-1, -1] = 1; T = np.diag(t); B0 = Dm + K@T + T@K - T@Dm@T
    out = f"N = {N}: background margin {eigh(((np.eye(N+1)+1j*T).T@H@(np.eye(N+1)+1j*T)).imag, np.diag(m), eigvals_only=True)[0]:+.4f};"
    for r in (1, 2, 3):
        V = eigh((B0+B0.T)/2, np.diag(m))[1][:, :r]; P = V.T*m[None, :]
        def f(a):
            C = np.eye(N+1) + 1j*(T + V@a.reshape(r, r)@P); Mx = (C.T@H@C).imag; return -eigh((Mx+Mx.T)/2, np.diag(m), eigvals_only=True)[0]
        res = minimize(f, np.zeros(r*r), method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-4, fatol=1e-7))
        out += f"  r = {r}: best margin {-res.fun:+.4f} (largest |A| entry {abs(res.x).max():.2f});"
    print(out)
