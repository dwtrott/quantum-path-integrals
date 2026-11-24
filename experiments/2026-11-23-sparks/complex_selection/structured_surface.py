"""Structured adapted surface for the index-one path, floating point:
z_j = x_j + i [ y_s,j + t_j xi_j + clip_B( (sigma - t_j) xi_j ) + a V_j clip_B2( P xi ) ],
uniform inner slope sigma, smooth outer tilt t(s) = 0.02 + 1.4 s (1 - s), and ONE smooth mode V (lowest eigenvector of
K + M relative to M, P = V^T M) with coefficient a.  sigma and a are fixed by maximising the local margin at N = 32.
Diagonal bent contours plus one clipped mode: every kick argument depends on its own coordinate and on one scalar."""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from nslice import Phi, grad, hess, setup
from diag_plus_low import saddle
def parts(N):
    z = saddle(N, -0.6, 0.3); d, w, K = setup(N); m = d*w; V = eigh(K + np.diag(m), np.diag(m))[1][:, 0]; V = V*np.sign(V.sum())
    s = np.linspace(0, 1, N+1); return z, m, V, 0.02 + 1.4*s*(1-s), hess(z, N)
def local(N, sigma, a):
    z, m, V, t, H = parts(N); C = np.eye(N+1) + 1j*(sigma*np.eye(N+1) + a*np.outer(V, V*m)); Mx = (C.T@H@C).imag
    return eigh((Mx+Mx.T)/2, np.diag(m), eigvals_only=True)[0]
res = minimize(lambda p: -local(32, np.clip(p[0], 0.02, 0.95), p[1]), [0.8, -1.5], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-7))
res2 = minimize(lambda p: -local(32, np.clip(p[0], 0.02, 0.95), p[1]), [0.8, 1.5], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-7))
if res2.fun < res.fun: res = res2
sigma, a = float(np.clip(res.x[0], 0.02, 0.95)), float(res.x[1]); print(f"sigma = {sigma:.4f}, a = {a:.4f} (fixed from N = 32)")
for N in (16, 32, 64, 128):
    z, m, V, t, H = parts(N); qb_pb = None; line = f"N = {N:3d}: local margin {local(N, sigma, a):+.4f};"
    if N <= 64:
        from nslice import classical
        x, p = classical(N); qb, pb = x[-1]-0.6, p+0.3; xs, ys, h0 = z.real, z.imag, Phi(z, N, qb, pb).imag; P = V*m
        for B, B2 in ((0.1, 0.1), (0.2, 0.2), (0.4, 0.4)):
            def zmap(x):
                xi = x-xs; l1 = (sigma-t)*xi; q = P@xi; return x + 1j*(ys + t*xi + np.clip(l1, -B, B) + a*V*np.clip(q, -B2, B2)), (abs(l1) < B), (abs(q) < B2)
            def h(x): return Phi(zmap(x)[0], N, qb, pb).imag - h0
            def dh(x):
                zz, m1, m2 = zmap(x); gr = grad(zz, N, qb, pb); dy = t + m1*(sigma-t); return gr.imag + dy*gr.real + (a*m2)*P*(V@gr.real)
            rng = np.random.default_rng(21); lo = np.inf; arg = xs
            for k in range(400 if N < 64 else 200):
                x0 = xs + rng.normal(0, rng.choice([0.2, 0.7, 2.0, 5.0]), N+1); r = minimize(h, x0, jac=dh, method='L-BFGS-B')
                if r.fun < lo: lo, arg = r.fun, r.x
            line += f"  B = {B}: min h {lo:+.5f} (dist {np.linalg.norm(arg-xs):.2f});"
    print(line)
