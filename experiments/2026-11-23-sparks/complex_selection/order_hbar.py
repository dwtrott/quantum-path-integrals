"""Toward a finite-hbar statement, floating point.
(c) Order-hbar coefficient a_1 of I/L = 1 + a_1 hbar + O(hbar^2), from the complex path's Green function G = H^(-1):
      a_1 = (i g / 8) sum_j m_j W''''(z_j) G_jj^2 + (i g^2 / 72) sum_jk m_j m_k W'''(z_j) W'''(z_k) (9 G_jj G_jk G_kk + 6 G_jk^3),
    compared with the FFT reference: (I_ref/L - 1)/hbar and the second-order residual (I_ref/L - 1 - a_1 hbar)/hbar^2.
(a) Ratio rho = sqrt(|det G_C| / det Re G_C), G_C = -i C^T H C: the factor lost by taking absolute values of the Gaussian on a plane C.
(b) Fredholm-type ratio sqrt(det Q / det(Q - c M)) on the steepest descent tangent plane, Q = R^T(-iH)R = I in R-coordinates,
    i.e. prod_k (1 - c/lambda_k)^(-1/2) with lambda_k the generalized eigenvalues of (Re-part quadratic form on the plane, M)."""
import numpy as np
from scipy.linalg import eigh, sqrtm, inv
from nslice import g, Phi, hess, setup, classical, newton, reference, W
W3 = lambda z: (3*z - z**3)*W(z); W4 = lambda z: (z**4 - 6*z*z + 3)*W(z)
sigma, a = 0.8291, -1.6247
def saddle(N, dq, dp):
    x, p = classical(N); z = x.astype(complex)
    for s in np.linspace(0, 1, 21)[1:]: z = newton(z, N, x[-1]+s*dq, p+s*dp, it=12)
    return newton(z, N, x[-1]+dq, p+dp), x[-1]+dq, p+dp
def plane(N, kind):
    s = np.linspace(0, 1, N+1); d, w, K = setup(N); m = d*w
    if kind == "smooth": return np.diag(1 + 1j*(0.02 + 1.4*s*(1-s)))
    if kind == "uniform": return (1+0.1j)*np.eye(N+1)
    return (1+1j*sigma)*np.eye(N+1) + (1j*a/2)*np.outer(np.ones(N+1), m)
def lead(z, N, qb, pb, hb, C):
    d = 2.0/N; G = -1j*C.T@hess(z, N)@C; ev = np.linalg.eigvals(G); sgn, ld = np.linalg.slogdet(C)
    return np.exp((N+1)/2*np.log(2*np.pi*hb) + np.log(sgn) + ld - 0.5*np.sum(np.log(ev)) + 1j*Phi(z, N, qb, pb)/hb - N/2*np.log(2j*np.pi*hb*d))
def a1(z, N):
    d, w, K = setup(N); m = d*w; G = inv(hess(z, N)); gd = np.diag(G); w3 = m*W3(z); w4 = m*W4(z)
    return 1j*g/8*np.sum(w4*gd**2) + 1j*g*g/72*(9*(w3*gd)@G@(w3*gd) + 6*w3@(G**3)@w3)
for label, dq, dp, kind in (("classical target", 0, 0, "smooth"), ("displaced (+0.4, +0.4)", 0.4, 0.4, "smooth"), ("displaced (-0.6, +0.3), index one", -0.6, 0.3, "one-mode")):
    print(f"\n{label}")
    print("    N | a_1 | (I_ref/L - 1)/hbar at 1/20, 1/40, 1/80 | |I_ref/L - 1 - a_1 hbar|/hbar^2 at 1/20, 1/40, 1/80")
    for N in (16, 32, 64, 128):
        z, qb, pb = saddle(N, dq, dp); C = plane(N, kind); c1 = a1(z, N); r1 = []; r2 = []
        for hb in (0.05, 0.025, 0.0125):
            r = reference(N, qb, pb, hb, M=2**15)/lead(z, N, qb, pb, hb, C) - 1; r1.append(f"{(r/hb).real:+.3f}{(r/hb).imag:+.3f}i"); r2.append(f"{abs(r - c1*hb)/hb**2:.2f}")
        print(f"  {N:4d} | {c1.real:+.4f}{c1.imag:+.4f}i | " + ", ".join(r1) + " | " + ", ".join(r2))
print("\n(a), (b) at the classical target")
print("    N | rho uniform tilt 0.1 | rho smooth tilt | rho steepest descent plane | Fredholm ratio c = 0.5 | c = 1.0 | smallest lambda")
for N in (16, 32, 64, 128, 256):
    z, qb, pb = saddle(N, 0, 0); H = hess(z, N); d, w, K = setup(N); m = d*w; out = []
    for kind in ("uniform", "smooth"):
        C = plane(N, kind); G = -1j*C.T@H@C; out.append(np.exp(0.5*(np.linalg.slogdet(G)[1] - np.linalg.slogdet((G+G.conj().T).real/2)[1])))
    R = inv(sqrtm(-1j*H)); G = -1j*R.T@H@R; out.append(np.exp(0.5*(np.linalg.slogdet(G)[1] - np.linalg.slogdet((G+G.conj().T).real/2)[1])))
    X = R.real; lam = 1/eigh(X.T@np.diag(m)@X, eigvals_only=True)      # on the plane z = z_s + R u, Gaussian form is |u|^2; ||Re z||_M^2 = u^T X^T M X u
    lam = np.sort(lam); fr = [np.prod(1/np.sqrt(1 - c/lam[lam > c])) if np.all(lam > c) else float('nan') for c in (0.5, 1.0)]
    print(f"  {N:4d} | {out[0]:.3e} | {out[1]:.3e} | {out[2]:.6f} | {fr[0]:.4f} | {fr[1]:.4f} | {lam[0]:.4f}")
