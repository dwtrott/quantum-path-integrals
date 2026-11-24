"""N-slice test of the tilted-plane rule (floating point).  Wavepacket scattering on a Gaussian barrier V = g exp(-x^2/2).
Phi(x_0..x_N) = i(x_0-q0)^2/2 + p0 x_0 + sum (x_{j+1}-x_j)^2/(2d) - d sum' V(x_j) + i(x_N-qb)^2/2 - pb x_N   (trapezoid weights).
With the factor (2 pi i hbar d)^(-N/2) the integral equals <coherent(qb,pb)| Strang split-step^N |packet(q0,p0)> exactly,
so the reference amplitude is computed by a one-dimensional FFT propagation with the same N.
For each N: saddle, leading term versus reference, and the selection margin
  mu_N = smallest generalized eigenvalue of  (1-t^2) D + 2 t K + diag(d w_j ell_j)  relative to diag(d w_j),
where K is the kinetic matrix, D = diag(1,0,..,0,1), T = t on every slice, and ell_j = inf_u R_j(u) is Nova's scalar remainder bound
(here estimated on a grid, not enclosed)."""
import numpy as np, sys
from scipy.linalg import eigh
g, q0, p0, Ttot = 3.0, -2.0, 2.0, 2.0
W = lambda z: np.exp(-z*z/2); W1 = lambda z: -z*W(z); W2 = lambda z: (z*z-1)*W(z)
def setup(N):
    d = Ttot/N; w = np.ones(N+1); w[0] = w[-1] = 0.5
    K = (np.diag(np.r_[1, 2*np.ones(N-1), 1]) - np.diag(np.ones(N), 1) - np.diag(np.ones(N), -1))/d
    return d, w, K
def Phi(z, N, qb, pb):
    d, w, K = setup(N)
    return 1j*(z[0]-q0)**2/2 + p0*z[0] + np.sum((z[1:]-z[:-1])**2)/(2*d) - d*g*np.sum(w*W(z)) + 1j*(z[-1]-qb)**2/2 - pb*z[-1]
def grad(z, N, qb, pb):
    d, w, K = setup(N); gr = K@z - d*g*w*W1(z); gr = gr.astype(complex); gr[0] += 1j*(z[0]-q0) + p0; gr[-1] += 1j*(z[-1]-qb) - pb; return gr
def hess(z, N):
    d, w, K = setup(N); H = K.astype(complex) - np.diag(d*g*w*W2(z)); H[0, 0] += 1j; H[-1, -1] += 1j; return H
def newton(z, N, qb, pb, it=60):
    for _ in range(it): z = z - np.linalg.solve(hess(z, N), grad(z, N, qb, pb))
    return z
def classical(N):
    d, w, K = setup(N); x = np.zeros(N+1); x[0] = q0; p = p0 + d*0.5*g*x[0]*W(x[0])          # force = -V' = g x W
    for j in range(N):
        x[j+1] = x[j] + d*p; p = p + d*w[j+1]*g*x[j+1]*W(x[j+1])
    return x, p
def reference(N, qb, pb, hbar, L=40.0, M=2**14):
    d = Ttot/N; x = (np.arange(M)-M/2)*(2*L/M); k = 2*np.pi*np.fft.fftfreq(M, 2*L/M)
    psi = np.exp(-(x-q0)**2/(2*hbar) + 1j*p0*x/hbar); half = np.exp(-1j*d*g*W(x)/(2*hbar)); drift = np.exp(-1j*hbar*k*k*d/2)
    for _ in range(N): psi = half*np.fft.ifft(drift*np.fft.fft(half*psi))
    return np.sum(np.exp(-(x-qb)**2/(2*hbar) - 1j*pb*x/hbar)*psi)*(2*L/M)
def lead(z, N, qb, pb, hbar):
    d = Ttot/N; ev = np.linalg.eigvals(-1j*hess(z, N))
    logL = (N+1)/2*np.log(2*np.pi*hbar) - 0.5*np.sum(np.log(ev)) + 1j*Phi(z, N, qb, pb)/hbar - N/2*np.log(2j*np.pi*hbar*d)
    return np.exp(logL)
us = np.r_[-np.geomspace(60, 1e-3, 4000), np.geomspace(1e-3, 60, 4000)]
def margin(z, N, t):
    d, w, K = setup(N); c = 1+1j*t; ell = np.zeros(N+1)
    for j in range(N+1):
        R = 2*np.imag(-g*(W(z[j]+c*us) - W(z[j]) - c*us*W1(z[j])))/us**2; ell[j] = R.min()
    D = np.zeros((N+1, N+1)); D[0, 0] = D[-1, -1] = 1
    Mx = (1-t*t)*D + 2*t*K + np.diag(d*w*ell)
    return eigh(Mx, np.diag(d*w), eigvals_only=True)[0], ell.min()
if __name__ == "__main__":
    for label, dq, dp in (("classical target", 0.0, 0.0), ("displaced target (+0.4, +0.4)", 0.4, 0.4), ("displaced target (-0.6, +0.3)", -0.6, 0.3)):
        print(f"\n{label}")
        print("   N |  qb       pb     | max|Im z| | Im Phi_s | |I_ref| hbar=.05 | rel err .05 | rel err .025 | best t | margin mu_N | min ell")
        for N in (4, 8, 16, 32, 64, 128):
            x, p = classical(N); qb, pb = x[-1] + dq, p + dp; z = x.astype(complex)
            for s in np.linspace(0, 1, 21)[1:]: z = newton(z, N, x[-1] + s*dq, p + s*dp, it=12)
            z = newton(z, N, qb, pb); errs = []
            for hb in (0.05, 0.025):
                I = reference(N, qb, pb, hb); Ls = lead(z, N, qb, pb, hb); errs.append(min(abs(I/Ls-1), abs(I/Ls+1)))
            best = max((margin(z, N, t) + (t,) for t in (0.02, 0.05, 0.1, 0.2, 0.3)), key=lambda r: r[0])
            print(f" {N:4d} | {qb:+.4f} {pb:+.4f} | {np.max(abs(z.imag)):.4f} | {Phi(z, N, qb, pb).imag:.5f} | {abs(reference(N, qb, pb, 0.05)):.4e} | {errs[0]:.4f} | {errs[1]:.4f} | {best[2]:.2f} | {best[0]:+.4f} | {best[1]:+.3f}")
