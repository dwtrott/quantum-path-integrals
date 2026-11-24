"""Exploration toward a selection rule (floating point).  Two-kick wavepacket model of correspondence 018.
The saddles depend on the target (b, pb) only through beta_bar = b - i pb, holomorphically.  So degenerate saddles
(turning points) are isolated points of the real target plane and the selection is a Stokes geometry on that plane.
1. Continue the real classical saddle from the classical target to every target on a grid.
2. Compare the exact amplitude with the leading term of the continued saddle.
3. Locate the turning points."""
import numpy as np, sys
a, p0, g = -2.0, 1.0, 3.0
def grad(z, b, pb):
    x1, x2 = z
    return np.array([(x1-a) - (1j-1)*(p0-(x1-a))/2 - (x2-x1) + g*x1*np.exp(-x1*x1/2),
                     (x2-x1) + (x2-b) + (1j-1)*(pb+(x2-b))/2 + g*x2*np.exp(-x2*x2/2)])
def hess(z):
    x1, x2 = z
    return np.array([[2 + (1j-1)/2 + g*(1-x1*x1)*np.exp(-x1*x1/2), -1], [-1, 2 + (1j-1)/2 + g*(1-x2*x2)*np.exp(-x2*x2/2)]])
def Phi(x1, x2, b, pb):
    return ((x1-a)**2/2 + (1j-1)*(p0-(x1-a))**2/4 + (x2-x1)**2/2 + (x2-b)**2/2 + (1j-1)*(pb+(x2-b))**2/4 - g*(np.exp(-x1*x1/2) + np.exp(-x2*x2/2)))
def newton(z, b, pb, it=40):
    for _ in range(it): z = z - np.linalg.solve(hess(z), grad(z, b, pb))
    return z
def exact(b, pb, hbar, half=7.0, n=1401):
    c1, c2 = a+p0, b-pb; w = half*np.sqrt(2*hbar)
    x1 = np.linspace(c1-w, c1+w, n); x2 = np.linspace(c2-w, c2+w, n); h = x1[1]-x1[0]
    X1, X2 = np.meshgrid(x1, x2, indexing='ij')
    return np.sum(np.exp(1j*Phi(X1, X2, b, pb)/hbar))*h*h
def lead(z, b, pb, hbar):
    M = -1j*hess(z); ev = np.linalg.eigvals(M)
    return 2*np.pi*hbar/np.prod(np.sqrt(ev))*np.exp(1j*Phi(z[0], z[1], b, pb)/hbar)
# classical target: real trajectory from the packet centre
x1c = a + p0; p1 = p0 + g*x1c*np.exp(-x1c*x1c/2); x2c = x1c + p1; p2 = p1 + g*x2c*np.exp(-x2c*x2c/2); bc, pbc = x2c + p2, p2
zc = newton(np.array([x1c, x2c], dtype=complex), bc, pbc)
print(f"classical target (b, pb) = ({bc:.4f}, {pbc:.4f}); real saddle ({zc[0]:.5f}, {zc[1]:.5f}), |grad| = {np.linalg.norm(grad(zc, bc, pbc)):.1e}")
def continue_to(b, pb, steps=400):
    z = zc.copy(); mind = np.inf
    for t in np.linspace(0, 1, steps+1)[1:]:
        z = newton(z, bc + t*(b-bc), pbc + t*(pb-pbc), it=8)
        mind = min(mind, abs(np.linalg.det(hess(z))))
    return newton(z, b, pb), mind
hbar = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
bs = np.linspace(-7, 5, 25); ps = np.linspace(-5, 4, 19)
print(f"\nhbar = {hbar}.  Entries: '.' relative deviation |I/L - 1| < 0.1, 'o' < 0.5, 'X' >= 0.5 (continued saddle does not give the amplitude), '-' amplitude below 1e-9 (not resolved)")
print("columns: b from -7 to 5; rows: pb from 4 (top) to -5")
rows = []; devs = {}
for pb in ps[::-1]:
    row = ""
    for b in bs:
        z, mind = continue_to(b, pb, steps=120); I = exact(b, pb, hbar); L = lead(z, b, pb, hbar)
        if abs(I) < 1e-9: row += "-"; continue
        d = min(abs(I/L - 1), abs(I/L + 1)); devs[(b, pb)] = (d, Phi(z[0], z[1], b, pb).imag, abs(I))
        row += "." if d < 0.1 else ("o" if d < 0.5 else "X")
    rows.append(row); print(f"pb={pb:+5.1f} " + row)
# 3. turning points: grad = 0 and det hess = 0 in (z1, z2, beta_bar)
def F(v):
    z = v[:2]; bb = v[2]; b, pb = bb.real, -bb.imag     # for complex beta_bar use the analytic form below
    x1, x2 = z
    g1 = (x1-a) - (1j-1)*(p0-(x1-a))/2 - (x2-x1) + g*x1*np.exp(-x1*x1/2)
    g2 = (x2-x1) + x2 + (1j-1)*x2/2 - (1+1j)/2*bb + g*x2*np.exp(-x2*x2/2)
    return np.array([g1, g2, np.linalg.det(hess(z))])
def numjac(v, eps=1e-7):
    J = np.zeros((3, 3), dtype=complex); f0 = F(v)
    for k in range(3):
        dv = np.zeros(3, dtype=complex); dv[k] = eps; J[:, k] = (F(v+dv) - F(v-dv))/(2*eps)
    return J
rng = np.random.default_rng(3); tps = []
for _ in range(4000):
    v = np.array([rng.uniform(-3, 3) + 1j*rng.uniform(-1.5, 1.5), rng.uniform(-3, 3) + 1j*rng.uniform(-1.5, 1.5), rng.uniform(-7, 5) - 1j*rng.uniform(-5, 4)])
    ok = True
    with np.errstate(all='ignore'):
        for _ in range(60):
            try: v = v - np.linalg.solve(numjac(v), F(v))
            except np.linalg.LinAlgError: ok = False; break
            if not np.all(np.isfinite(v)) or np.max(abs(v)) > 40: ok = False; break
        if ok and np.linalg.norm(F(v)) < 1e-9 and all(abs(v[2]-w[2]) > 1e-5 or np.linalg.norm(v[:2]-w[:2]) > 1e-5 for w in tps): tps.append(v)
tps = [v for v in tps if -7.5 <= v[2].real <= 5.5 and -4.5 <= v[2].imag <= 5.5 and max(abs(v[:2].imag)) < 1.2]
print(f"\nturning points (degenerate saddles) with target in the mapped window and |Im z| < 1.2: {len(tps)}")
for v in sorted(tps, key=lambda v: abs(v[2] - (bc - 1j*pbc)))[:14]:
    b, pb = v[2].real, -v[2].imag
    print(f"   target (b, pb) = ({b:+.3f}, {pb:+.3f})  distance from classical target {abs(v[2]-(bc-1j*pbc)):.3f}  saddle ({v[0]:.3f}, {v[1]:.3f})  Im Phi = {Phi(v[0], v[1], b, pb).imag:+.3f}")
