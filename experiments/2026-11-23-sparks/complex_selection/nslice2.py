"""Slice-dependent tilt on the N-slice barrier model (floating point).  Tilt t_j = t(s_j), t(s) quadratic in s plus separate
endpoint values.  Margin = smallest generalized eigenvalue of  D + K T + T K - T D T + diag(d w_j ell_j)  relative to diag(d w_j),
ell_j = inf_u R_j(u) with c_j = 1 + i t_j (grid estimate, not an enclosure).
Leading term with the oriented Gaussian branch on the tilted plane: det(C) / sqrt(det(-i C H C)), principal roots of the
eigenvalues, valid when the real part of -i C H C is positive definite.  No sign fitting."""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from nslice import g, q0, p0, Ttot, W, W1, setup, Phi, hess, newton, classical, reference, us
def tilt(par, N):
    a, b, c, e0, eN = par; s = np.linspace(0, 1, N+1); t = a + b*s + c*s*s; t[0], t[-1] = e0, eN
    return np.clip(t, -0.9, 0.9)
def margin(z, N, t):
    d, w, K = setup(N); ell = np.zeros(N+1)
    for j in range(N+1):
        c = 1+1j*t[j]; ell[j] = (2*np.imag(-g*(W(z[j]+c*us) - W(z[j]) - c*us*W1(z[j])))/us**2).min()
    D = np.zeros((N+1, N+1)); D[0, 0] = D[-1, -1] = 1; T = np.diag(t)
    return eigh(D + K@T + T@K - T@D@T + np.diag(d*w*ell), np.diag(d*w), eigvals_only=True)[0]
def lead_oriented(z, N, qb, pb, hbar, t):
    d = Ttot/N; C = np.diag(1+1j*t); Mx = -1j*C@hess(z, N)@C
    pd = np.linalg.eigvalsh((Mx+Mx.conj().T)/2)[0] > 0; ev = np.linalg.eigvals(Mx)
    logL = (N+1)/2*np.log(2*np.pi*hbar) + np.sum(np.log(1+1j*t)) - 0.5*np.sum(np.log(ev)) + 1j*Phi(z, N, qb, pb)/hbar - N/2*np.log(2j*np.pi*hbar*d)
    return np.exp(logL), pd
def saddle(N, dq, dp):
    x, p = classical(N); z = x.astype(complex)
    for s in np.linspace(0, 1, 21)[1:]: z = newton(z, N, x[-1]+s*dq, p+s*dp, it=12)
    return newton(z, N, x[-1]+dq, p+dp), x[-1]+dq, p+dp
for label, dq, dp in (("classical target", 0, 0), ("displaced (+0.4, +0.4)", 0.4, 0.4), ("displaced (-0.6, +0.3)", -0.6, 0.3)):
    z, qb, pb = saddle(32, dq, dp); best = None
    for start in ([0.2, 0, 0, 0.2, 0.2], [0.1, 0.3, 0, 0.0, 0.3], [0.4, -0.3, 0, 0.3, 0.0], [0.3, 0, 0, 0, 0], [0.5, 0, -0.4, 0.2, 0.2]):
        r = minimize(lambda par: -margin(z, 32, tilt(par, 32)), start, method='Nelder-Mead', options=dict(maxiter=400, xatol=1e-3, fatol=1e-5))
        if best is None or r.fun < best.fun: best = r
    print(f"\n{label}: optimised tilt parameters (a, b, c, t_0, t_N) = {np.round(best.x, 3)}; im z profile at N=32: first {z[0].imag:+.3f}, max |Im z| {np.max(abs(z.imag)):.3f} at slice {int(np.argmax(abs(z.imag)))}, last {z[-1].imag:+.3f}")
    print("   N | margin | oriented branch valid | signed rel err hbar=.05 | hbar=.025")
    for N in (8, 16, 32, 64, 128):
        z, qb, pb = saddle(N, dq, dp); t = tilt(best.x, N); errs = []
        for hb in (0.05, 0.025):
            L, pd = lead_oriented(z, N, qb, pb, hb, t); errs.append(abs(reference(N, qb, pb, hb)/L - 1))
        print(f" {N:4d} | {margin(z, N, t):+.4f} | {pd} | {errs[0]:.4f} | {errs[1]:.4f}")
