"""Oriented leading term for the index-one path (third target) on the one-mode plane, floating point.
Local plane w = u v + b mu, u = 1 + i sigma, b = i a, mu = m^T v / 2, i.e. C = u I + (b/2) 1 m^T  (Nova 039: the smoothest mode
is the constant vector).  On this plane Im(C^T H C) is positive definite, so the Gaussian integral is
det(C) (2 pi hbar)^((N+1)/2) / sqrt(det(-i C^T H C)) with principal square roots of the eigenvalues: no sign is fitted.
Reference: FFT Strang propagation with the same N (periodic box approximation of the exact identity).
Also reports the mass-normalised local margin and its Richardson extrapolation in N."""
import numpy as np
from scipy.linalg import eigh
from nslice import Phi, hess, setup, classical, newton, reference, Ttot
sigma, a = 0.8291, -1.6247
def saddle(N, dq, dp):
    x, p = classical(N); z = x.astype(complex)
    for s in np.linspace(0, 1, 21)[1:]: z = newton(z, N, x[-1]+s*dq, p+s*dp, it=12)
    return newton(z, N, x[-1]+dq, p+dp), x[-1]+dq, p+dp
marg = {}
print("   N | local margin | Re(-i C^T H C) > 0 | I_ref/L - 1 at hbar = 1/20 | 1/40 | 1/80")
for N in (16, 32, 64, 128, 256):
    z, qb, pb = saddle(N, -0.6, 0.3); d, w, K = setup(N); m = d*w; H = hess(z, N)
    C = (1+1j*sigma)*np.eye(N+1) + (1j*a/2)*np.outer(np.ones(N+1), m); G = -1j*C.T@H@C
    herm = (G+G.conj().T)/2; marg[N] = eigh(herm.real, np.diag(m), eigvals_only=True)[0]; pd = np.linalg.eigvalsh(herm)[0] > 0
    ev = np.linalg.eigvals(G); sgn, logdetC = np.linalg.slogdet(C); out = []
    for hb in (0.05, 0.025, 0.0125):
        logL = (N+1)/2*np.log(2*np.pi*hb) + np.log(sgn) + logdetC - 0.5*np.sum(np.log(ev)) + 1j*Phi(z, N, qb, pb)/hb - N/2*np.log(2j*np.pi*hb*d)
        r = reference(N, qb, pb, hb, M=2**15)/np.exp(logL) - 1; out.append(f"{r.real:+.4f}{r.imag:+.4f}i")
    print(f" {N:4d} | {marg[N]:+.5f} | {pd} | " + " | ".join(out))
print(f"Richardson extrapolation of the margin from N = 128, 256: {(4*marg[256]-marg[128])/3:.5f};  successive differences {marg[32]-marg[16]:.5f}, {marg[64]-marg[32]:.5f}, {marg[128]-marg[64]:.5f}, {marg[256]-marg[128]:.5f}")
