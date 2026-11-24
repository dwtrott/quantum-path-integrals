"""Third target (Morse index one), floating point.  All margins are normalised by the mass matrix diag(d w_j).
(1) Diagonal tilts: best local margin over ALL per-slice tilts, and over smooth positive tilts (degree-4 polynomial in the
    time fraction, 0.02 <= t <= 0.9).  Corrects the unnormalised numbers of correspondence 035, Section 4.
(2) Clipped-thimble surface with the smooth tilt 0.02 + 1.4 s (1 - s): local margin and multistart global minimum of h."""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from nslice import hess, setup
from nslice2 import saddle
from clipped_thimble import run, surface
def local_margin(H, t, N):
    d, w, K = setup(N); c = 1+1j*t; return eigh((c[:, None]*H*c[None, :]).imag, np.diag(d*w), eigvals_only=True)[0]
def poly(par, N):
    s = np.linspace(0, 1, N+1); return np.clip(sum(p*s**k for k, p in enumerate(par)), 0.02, 0.9)
rng = np.random.default_rng(4)
print("(1) diagonal tilts, local margin at the saddle")
zb, _, _ = saddle(32, -0.6, 0.3); Hb = hess(zb, 32); bestp = None
for k in range(10):
    p0 = np.r_[rng.uniform(0.05, 0.5), rng.uniform(-1.5, 1.5, 4)]
    r = minimize(lambda p: -local_margin(Hb, poly(p, 32), 32), p0, method='Nelder-Mead', options=dict(maxiter=3000, xatol=1e-4, fatol=1e-7))
    if bestp is None or r.fun < bestp.fun: bestp = r
print(f"   smooth positive tilt optimised at N = 32: coefficients {np.round(bestp.x, 3)}")
for N in (8, 16, 32, 64, 128):
    z, _, _ = saddle(N, -0.6, 0.3); H = hess(z, N); best = np.inf
    if N <= 32:
        for k in range(8):
            t0 = rng.uniform(-0.3, 0.6, N+1) if k else np.full(N+1, 0.2)
            r = minimize(lambda t: -local_margin(H, np.clip(t, -0.95, 0.95), N), t0, method='Powell', options=dict(maxiter=20000, xtol=1e-4, ftol=1e-7)); best = min(best, r.fun)
    print(f"   N = {N:3d}: free per-slice tilt {(-best if N <= 32 else float('nan')):+.4f};  smooth positive tilt {local_margin(H, poly(bestp.x, N), N):+.4f};  tilt 0.02+1.4s(1-s) {local_margin(H, poly([0.02, 1.4, -1.4], N), N):+.4f}")
print("(2) clipped-thimble surface, tilt 0.02 + 1.4 s (1 - s)")
for N in (16, 32, 64):
    s = np.linspace(0, 1, N+1); t = 0.02 + 1.4*s*(1-s); z, _, _ = saddle(N, -0.6, 0.3); d, w, K = setup(N)
    zmap, G, S = surface(z, N, t, 0.1); C = np.eye(N+1) + 1j*S; M = (C.T@hess(z, N)@C).imag
    print(f"   N = {N}: mass-normalised local margin on the thimble tangent plane {eigh((M+M.T)/2, np.diag(d*w), eigvals_only=True)[0]:+.4f}")
    for B in (0.05, 0.1, 0.15, 0.2): run(N, -0.6, 0.3, "smooth", t, B, starts=600 if N < 64 else 250, seed=7)
