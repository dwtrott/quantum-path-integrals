"""Local version of hypothesis (iii): largest ratio (|u|^2 - 2 h(u)) / ||X u||_M^2 found over the box |u_k| <= beta on the
steepest descent tangent plane, classical target.  Multistart with bounds; floating point, not a certificate."""
import numpy as np
from scipy.linalg import sqrtm, inv, eigh
from scipy.optimize import minimize
from nslice import Phi, hess, setup
from order_hbar import saddle
rng = np.random.default_rng(5)
for N in (16, 32, 64):
    z, qb, pb = saddle(N, 0, 0); H = hess(z, N); d, w, K = setup(N); m = d*w; R = inv(sqrtm(-1j*H)); X = R.real; h0 = Phi(z, N, qb, pb).imag
    lam1 = 1/eigh(X.T@np.diag(m)@X, eigvals_only=True)[-1]; line = f"N = {N}: lambda_1 = {lam1:.4f};"
    def negratio(u):
        hv = Phi(z + R@u, N, qb, pb).imag - h0; return -(u@u - 2*hv)/((X@u)@(m*(X@u)))
    for beta in (0.05, 0.1, 0.2, 0.4, 0.8):
        best = -np.inf; arg = None
        for k in range(120 if N < 64 else 50):
            u0 = rng.uniform(-beta, beta, N+1); r = minimize(negratio, u0, method='L-BFGS-B', bounds=[(-beta, beta)]*(N+1))
            if -r.fun > best: best, arg = -r.fun, r.x
        line += f"  beta = {beta}: c = {best:.3f} (max|delta| {np.max(abs(R@arg)):.2f});"
    print(line)
