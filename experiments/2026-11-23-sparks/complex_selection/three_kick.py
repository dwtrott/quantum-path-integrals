"""Three-kick wavepacket model with a free interior variable x2 (floating point).
Phi = (x1-a)^2/2 + (i-1)(p0-(x1-a))^2/4 + (x2-x1)^2/2 + (x3-x2)^2/2 + (x3-b)^2/2 + (i-1)(pb+(x3-b))^2/4 - g sum_j exp(-x_j^2/2).
Checks: (1) on any flat plane Im Phi has no coercivity in x2, so the flat rule is void; (2) a tilted plane
z_j = z_s,j + (1 + i t_j) u_j selects the contributing saddle.  The exact amplitude is computed with x2 on the fixed
contour x2 = (1 + 0.3 i) u, which is a legitimate deformation and independent of any saddle."""
import numpy as np, itertools, sys
a, p0, g = -2.0, 1.0, 3.0
W = lambda z: np.exp(-z*z/2)
def Phi(z1, z2, z3, b, pb):
    return ((z1-a)**2/2 + (1j-1)*(p0-(z1-a))**2/4 + (z2-z1)**2/2 + (z3-z2)**2/2 + (z3-b)**2/2 + (1j-1)*(pb+(z3-b))**2/4 - g*(W(z1)+W(z2)+W(z3)))
def grad(z, b, pb):
    z1, z2, z3 = z
    return np.array([(z1-a) - (1j-1)*(p0-(z1-a))/2 - (z2-z1) + g*z1*W(z1), (z2-z1) - (z3-z2) + g*z2*W(z2),
                     (z3-z2) + (z3-b) + (1j-1)*(pb+(z3-b))/2 + g*z3*W(z3)])
def hess(z):
    k = [g*(1-w*w)*W(w) for w in z]
    return np.array([[2 + (1j-1)/2 + k[0], -1, 0], [-1, 2 + k[1], -1], [0, -1, 2 + (1j-1)/2 + k[2]]])
def lead(z, b, pb, hbar):
    ev = np.linalg.eigvals(-1j*hess(z)); return (2*np.pi*hbar)**1.5/np.prod(np.sqrt(ev))*np.exp(1j*Phi(*z, b, pb)/hbar)
def exact(b, pb, hbar, n=261):
    w = 7*np.sqrt(2*hbar); x1 = np.linspace(a+p0-w, a+p0+w, n); x3 = np.linspace(b-pb-w, b-pb+w, n); u = np.linspace(-9, 9, 1401); c = 1+0.3j
    X1, X3 = np.meshgrid(x1, x3, indexing='ij'); tot = 0
    for uu in u: tot += np.sum(np.exp(1j*Phi(X1, c*uu, X3, b, pb)/hbar))
    return tot*c*(x1[1]-x1[0])*(x3[1]-x3[0])*(u[1]-u[0])
rng = np.random.default_rng(5); starts = [rng.uniform(-4, 4, 3) + 1j*rng.uniform(-1.6, 1.6, 3) for _ in range(2500)]
def saddles(b, pb):
    out = []
    with np.errstate(all='ignore'):
        for z in starts:
            z = z.copy(); ok = True
            for _ in range(50):
                try: z = z - np.linalg.solve(hess(z), grad(z, b, pb))
                except np.linalg.LinAlgError: ok = False; break
                if not np.all(np.isfinite(z)) or np.max(abs(z)) > 30: ok = False; break
            if ok and np.linalg.norm(grad(z, b, pb)) < 1e-10 and np.max(abs(z.imag)) < 3 and all(np.linalg.norm(z-w) > 1e-6 for w in out): out.append(z)
    return out
us = np.linspace(-10, 10, 101); U = np.meshgrid(us, us, us, indexing='ij')
def gap(z, b, pb, t):
    with np.errstate(all='ignore'):
        h = Phi(z[0]+(1+1j*t[0])*U[0], z[1]+(1+1j*t[1])*U[1], z[2]+(1+1j*t[2])*U[2], b, pb).imag
    return Phi(*z, b, pb).imag - np.nanmin(h)
# classical target
x1 = a+p0; p = p0 + g*x1*W(x1); x2 = x1+p; p = p + g*x2*W(x2); x3 = x2+p; p = p + g*x3*W(x3); bc, pbc = x3+p, p
print(f"classical target (b, pb) = ({bc:.4f}, {pbc:.4f})")
targets = [(bc, pbc), (bc+1.0, pbc), (bc, pbc+1.0), (0.0, 0.0), (2.0, 1.0), (-2.0, -1.0)]
tgrid = [(t1, t2, t3) for t1 in (-0.1, 0, 0.1) for t2 in (0.05, 0.1, 0.2, 0.35) for t3 in (-0.1, 0, 0.1)]
for b, pb in targets:
    I1, I2 = exact(b, pb, 0.1), exact(b, pb, 0.05)
    S = [z for z in saddles(b, pb) if Phi(*z, b, pb).imag > -1e-6]; S.sort(key=lambda z: Phi(*z, b, pb).imag)
    with np.errstate(all='ignore'): fits = [min(abs(I2/lead(z, b, pb, 0.05)-1), abs(I2/lead(z, b, pb, 0.05)+1)) for z in S]
    k = int(np.argmin(fits)); z = S[k]; rate = (np.log(abs(I1)/0.1**1.5) - np.log(abs(I2)/0.05**1.5))/10
    gflat = gap(z, b, pb, (0, 0, 0)); G = [(gap(z, b, pb, t), t) for t in tgrid]; gb, tb = min(G)
    lower = [w for j, w in enumerate(S) if j < k]; fs = sum(1 for w in lower if min(gap(w, b, pb, t) for t in tgrid[::4]) < 2e-3)
    print(f"target ({b:+.3f},{pb:+.3f}): |I(0.05)| = {abs(I2):.3e}, saddles with Im Phi >= 0: {len(S)}, contributing saddle rank {k} by Im Phi, Im Phi_s = {Phi(*z, b, pb).imag:.4f}, measured rate {rate:.4f}, fit {fits[k]:.3f}")
    print(f"     flat plane: depth below saddle on the box = {gflat:.2f} (unbounded as the box grows);  best tilt {tb}: gap {gb:.5f} -> {'tilt-selected' if gb < 2e-3 else 'not selected on this tilt grid'};  lower saddles falsely selected: {fs} of {len(lower)}")
