"""Clipped-thimble surface for the Morse-index-one path (third target of 033), floating point.
Surface as a graph over real x:   z = x + i [ y_s + t o (x - x_s) + clip_B( (S - diag t)(x - x_s) ) ],
S = Im(R) Re(R)^(-1), R = (-i H)^(-1/2), so that near the saddle the surface is the tangent plane of the steepest descent
manifold and far away it is the diagonal tilted plane shifted by a bounded amount (clip componentwise to [-B, B]).
Question: is the saddle the global minimum of h(x) = Im Phi(z(x)) - Im Phi(z_s) on this surface?  Multistart minimisation."""
import numpy as np, sys
from scipy.linalg import sqrtm, inv
from scipy.optimize import minimize
from nslice import Phi, grad, hess
from nslice2 import saddle
def surface(z, N, t, B):
    H = hess(z, N); R = inv(sqrtm(-1j*H)); S = R.imag@inv(R.real); G = S - np.diag(t); xs, ys = z.real, z.imag
    def zmap(x):
        lin = G@(x-xs); return x + 1j*(ys + t*(x-xs) + np.clip(lin, -B, B)), (abs(lin) < B)
    return zmap, G, S
def run(N, dq, dp, tname, t, B, starts=500, seed=0):
    z, qb, pb = saddle(N, dq, dp); zmap, G, S = surface(z, N, t, B); h0 = Phi(z, N, qb, pb).imag; xs = z.real
    def h(x): return Phi(zmap(x)[0], N, qb, pb).imag - h0
    def dh(x):
        zz, mask = zmap(x); J = np.diag(t) + mask[:, None]*G; gr = grad(zz, N, qb, pb); return gr.imag + J.T@gr.real
    rng = np.random.default_rng(seed); lo = np.inf; arg = None
    for k in range(starts):
        x0 = xs + rng.normal(0, rng.choice([0.2, 0.7, 2.0, 5.0]), N+1); r = minimize(h, x0, jac=dh, method='L-BFGS-B')
        if r.fun < lo: lo, arg = r.fun, r.x
    # local Hessian of h at the saddle on the surface: Im((I+iS)^T H (I+iS))
    C = np.eye(N+1) + 1j*S; loc = np.linalg.eigvalsh(((C.T@hess(z, N)@C).imag + (C.T@hess(z, N)@C).imag.T)/2)[0]
    print(f"   N = {N:3d}, tilt {tname}, B = {B}: smallest h found {lo:+.6f} at distance {np.linalg.norm(arg-xs):.3f};  local Hessian min eig {loc:+.4f};  max |S - diag t| entry {abs(G).max():.2f}")
    return lo
if __name__ == "__main__":
    for N in (8, 16, 32):
        s = np.linspace(0, 1, N+1)
        for tname, t in (("0.2", np.full(N+1, 0.2)), ("0.02+1.4s(1-s)", 0.02 + 1.4*s*(1-s))):
            for B in (0.1, 0.3, 0.6, 1.0):
                run(N, -0.6, 0.3, tname, t, B, starts=300 if N < 32 else 150)
