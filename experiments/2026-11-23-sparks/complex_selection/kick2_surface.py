"""Float-level design test of a two-variable selection certificate (not a proof).
Surface: Gamma = { x + i (y_s + clip_B(K (x - x_s))) : x in R^2 }, K = Im(R) Re(R)^{-1}, R = (-i Phi''(z_s))^{-1/2},
so that Gamma is the exact tangent plane of the thimble near the saddle and stays in a strip of half-width ~B.
Local lemma on the square |u|_inf <= rho in the coordinates z = z_s + R u; tail by direct numerical integration."""
import numpy as np, sys
from scipy.linalg import sqrtm, inv
from scipy.special import erfc
from kick2_explore import make, saddles, lead, exact, g
def W3(z): return np.exp(-z*z/2)*(3*z - z**3)          # W''' for W = exp(-z^2/2)
def W4(z): return np.exp(-z*z/2)*(z**4 - 6*z*z + 3)
def design(b, pb, idx, Bclip, kappa, hbars, competitor=None, ngrid=1601, half=6.0):
    Phi, grad, hess = make(b, pb); S = saddles(grad, hess, Phi)
    pos = [z for z in S if Phi(z[0], z[1]).imag >= -1e-9]; zs = pos[idx]
    H = hess(zs); R = inv(sqrtm(-1j*H)); K = R.imag @ inv(R.real)
    xs, ys = zs.real, zs.imag
    assert np.linalg.det(R.real) > 0
    rows = np.linalg.norm(R, axis=1)
    # largest admissible rho: sqrt(2) rho M3 <= 3 kappa, and the square inside the unclipped region
    rho = 0.0
    for cand in np.linspace(0.02, 2.0, 100):
        uu = np.linspace(-cand, cand, 41); U1, U2 = np.meshgrid(uu, uu)
        Z = zs[:, None, None] + np.einsum('ij,jkl->ikl', R, np.array([U1, U2]))
        M3 = g*np.max(abs(W3(Z)))*np.sum(rows**3)
        yoff = np.einsum('ij,jkl->ikl', R.imag, np.array([U1, U2]))
        if np.sqrt(2)*cand*M3 <= 3*kappa and np.max(abs(yoff)) <= Bclip: rho = cand; M3k = M3; M4k = g*np.max(abs(W4(Z)))*np.sum(rows**4)
        else: break
    print(f" saddle z_s = ({zs[0]:.5f}, {zs[1]:.5f}), Im Phi_s = {Phi(zs[0],zs[1]).imag:.5f}; rho = {rho:.3f}, M3 = {M3k:.3f}, M4 = {M4k:.3f}")
    def surf(X1, X2):
        V = np.einsum('ij,jkl->ikl', K, np.array([X1 - xs[0], X2 - xs[1]]))
        Y = np.clip(V, -Bclip, Bclip)
        # Jacobian of y wrt x: rows of K kept where unclipped
        m1 = (abs(V[0]) < Bclip).astype(float); m2 = (abs(V[1]) < Bclip).astype(float)
        J = (1 + 1j*m1*K[0, 0])*(1 + 1j*m2*K[1, 1]) - (1j*m1*K[0, 1])*(1j*m2*K[1, 0])
        return X1 + 1j*(ys[0] + Y[0]), X2 + 1j*(ys[1] + Y[1]), J
    Rinv = inv(R.real)
    for hbar in hbars:
        w = half*np.sqrt(2*hbar) + 1.0
        x1 = np.linspace(xs[0]-w, xs[0]+w, ngrid); x2 = np.linspace(xs[1]-w, xs[1]+w, ngrid); h = x1[1]-x1[0]
        X1, X2 = np.meshgrid(x1, x2, indexing='ij'); Z1, Z2, J = surf(X1, X2)
        F = np.exp(1j*Phi(Z1, Z2)/hbar)*J
        Isurf = np.sum(F)*h*h; Iex = exact(Phi, b, pb, hbar)
        U = np.einsum('ij,jkl->ikl', Rinv, np.array([X1 - xs[0], X2 - xs[1]]))
        inside = (abs(U[0]) <= rho) & (abs(U[1]) <= rho)
        T = np.sum(abs(F[~inside]))*h*h
        ImP = Phi(Z1, Z2).imag; mu = ImP[~inside].min()
        L0 = 2*np.pi*hbar*np.linalg.det(R)*np.exp(1j*Phi(zs[0], zs[1])/hbar)
        r = 2*erfc(rho/np.sqrt(2*hbar)) + hbar*(M4k/3 + (2/3)*M3k**2/(1-kappa)**4)
        Bnd = abs(L0)*r + T
        line = f"  hbar={hbar}: |I|={abs(Iex):.4e}  surface-vs-real-plane check {abs(Isurf-Iex):.1e}  true error |I-L0|={abs(Iex-L0):.3e}  bound={Bnd:.3e} (local {abs(L0)*r:.3e}, tail {T:.3e}; mu-ImPhi_s={mu-Phi(zs[0],zs[1]).imag:.4f})  relative bound {Bnd/abs(L0):.3f}"
        if competitor is not None:
            zc = pos[competitor]; Lc = abs(lead(zc, Phi, hess, hbar)); line += f"  |competitor term|={Lc:.3e}  excluded: {Bnd < Lc}"
        print(line)
for (b, pb, idx, comp, label) in ((-3.68, -1.86, 0, 1, "allowed target; competitor saddle has Im Phi = 0.165"), (2.0, 1.0, 0, None, "forbidden target; leading saddle complex, Im Phi = 0.258")):
    for Bclip in (0.6, 1.0):
        print(f"\n== {label}; clip B = {Bclip}, kappa = 0.25")
        design(b, pb, idx, Bclip, 0.25, (0.05, 0.025, 0.0125, 0.00625), competitor=comp)
