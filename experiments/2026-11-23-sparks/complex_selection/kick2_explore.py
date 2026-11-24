"""Two kicks, wavepacket in and wavepacket out: a two-variable integral with complex coefficients.
I = int int exp(i Phi(x1,x2)/hbar) dx1 dx2,
Phi = phi_in(x1) + (x2-x1)^2/2 + phi_out(x2) - g (W(x1) + W(x2)),  W(x) = exp(-x^2/2),
phi_in(x) = (x-a)^2/2 + (i-1)(p0-(x-a))^2/4,   phi_out(x) = (x-b)^2/2 + (i-1)(pb+(x-b))^2/4.
On R^2, Im Phi = [(x1-a-p0)^2 + (x2-b+pb)^2]/4, so the integral converges absolutely.
Floating point exploration: exact value by the trapezoid rule; complex saddles by Newton; which subset reproduces I?"""
import numpy as np, itertools, sys
a, p0, g = -2.0, 1.0, 3.0
def make(b, pb):
    def Phi(x1, x2):
        return ((x1-a)**2/2 + (1j-1)*(p0-(x1-a))**2/4 + (x2-x1)**2/2 + (x2-b)**2/2 + (1j-1)*(pb+(x2-b))**2/4
                - g*(np.exp(-x1*x1/2) + np.exp(-x2*x2/2)))
    def grad(z):
        x1, x2 = z
        return np.array([(x1-a) - (1j-1)*(p0-(x1-a))/2 - (x2-x1) + g*x1*np.exp(-x1*x1/2),
                         (x2-x1) + (x2-b) + (1j-1)*(pb+(x2-b))/2 + g*x2*np.exp(-x2*x2/2)])
    def hess(z):
        x1, x2 = z
        return np.array([[1 + (1j-1)/2 + 1 + g*(1-x1*x1)*np.exp(-x1*x1/2), -1],
                         [-1, 1 + 1 + (1j-1)/2 + g*(1-x2*x2)*np.exp(-x2*x2/2)]])
    return Phi, grad, hess
def exact(Phi, b, pb, hbar, half=7.0, n=3201):
    c1, c2 = a+p0, b-pb; w = half*np.sqrt(2*hbar)
    x1 = np.linspace(c1-w, c1+w, n); x2 = np.linspace(c2-w, c2+w, n); h = x1[1]-x1[0]
    X1, X2 = np.meshgrid(x1, x2, indexing='ij')
    return np.sum(np.exp(1j*Phi(X1, X2)/hbar))*h*h
def saddles(grad, hess, Phi, tries=6000, seed=0):
    rng = np.random.default_rng(seed); found = []
    for _ in range(tries):
        z = rng.uniform(-3.5, 3.5, 2) + 1j*rng.uniform(-2.5, 2.5, 2)
        ok = True
        with np.errstate(all='ignore'):
            for _ in range(80):
                try: dz = np.linalg.solve(hess(z), grad(z))
                except np.linalg.LinAlgError: ok = False; break
                z = z - dz
                if not np.all(np.isfinite(z)) or np.max(abs(z)) > 30: ok = False; break
            if ok and np.linalg.norm(grad(z)) < 1e-10 and np.max(abs(z.real)) < 5 and np.max(abs(z.imag)) < 4 and all(np.linalg.norm(z-w) > 1e-6 for w in found):
                found.append(z)
    return sorted(found, key=lambda z: Phi(z[0], z[1]).imag)
def lead(z, Phi, hess, hbar):
    M = -1j*hess(z); ev = np.linalg.eigvals(M)
    return 2*np.pi*hbar/np.prod(np.sqrt(ev+0j))*np.exp(1j*Phi(z[0], z[1])/hbar)   # sign fixed below by comparison
if __name__ == "__main__":
    for (b, pb, label) in ((-3.68, -1.86, "classically allowed target"), (2.0, 1.0, "transmitted target, classically forbidden"), (0.0, 0.0, "intermediate target")):
        Phi, grad, hess = make(b, pb); S = saddles(grad, hess, Phi)
        print(f"\n=== target (b, pb) = ({b}, {pb}): {label};  {len(S)} saddles found in the search box")
        neg = [z for z in S if Phi(z[0], z[1]).imag < -1e-9]; pos = [z for z in S if Phi(z[0], z[1]).imag >= -1e-9]
        print(f"  saddles with Im Phi < 0: {len(neg)} (most negative Im Phi = {min([Phi(z[0],z[1]).imag for z in neg], default=0):.2f});  with Im Phi >= 0: {len(pos)}")
        for z in pos[:5]: print(f"    z = ({z[0]:.4f}, {z[1]:.4f})   Im Phi = {Phi(z[0],z[1]).imag:.5f}")
        for hbar in (0.2, 0.1, 0.07):
            I = exact(Phi, b, pb, hbar); I2 = exact(Phi, b, pb, hbar, n=2401)
            cands = pos[:4]; best = None
            for k in range(1, len(cands)+1):
                for sub in itertools.combinations(range(len(cands)), k):
                    for signs in itertools.product((1, -1), repeat=k):
                        err = abs(I - sum(s*lead(cands[i], Phi, hess, hbar) for s, i in zip(signs, sub)))
                        if best is None or err < best[0]*0.999: best = (err, sub, signs)
            s1 = min(abs(I - lead(pos[0], Phi, hess, hbar)), abs(I + lead(pos[0], Phi, hess, hbar)))
            print(f"  hbar={hbar}: |I| = {abs(I):.4e} (grid check {abs(I-I2):.1e});  lowest saddle alone: error {s1:.3e} (relative {s1/abs(I):.3e});  best subset {best[1]} signs {best[2]}: error {best[0]:.3e}")
