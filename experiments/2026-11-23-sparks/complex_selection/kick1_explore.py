"""One kick, wavepacket in, position out.  Units m = d = 1, packet width sigma0 = 1 (in units of sqrt(hbar)).
Amplitude(b) = int exp(i phi(x)/hbar) dx,  phi(x) = (x-a)^2/2 + (i-1)(p0-(x-a))^2/4 + (b-x)^2/2 - g exp(-x^2/2).
All coefficients complex; on the real axis Im phi = (x-a-p0)^2/4 >= 0, so the integral converges absolutely.
Exploration in floating point / mpmath: which complex saddles reproduce the exact amplitude?"""
import numpy as np, mpmath as mp, itertools, sys
a, p0, g = -2.0, 1.0, 3.0
def coeffs(b):
    lam = 1.5 + 0.5j
    # phi = lam x^2/2 - c x + k0 - g exp(-x^2/2); expand the quadratic part
    # (x-a)^2/2 + (i-1)(p0+a-x)^2/4 + (b-x)^2/2
    q = p0 + a
    c = a + (1j-1)*q/2 + b              # minus the linear coefficient
    k0 = a*a/2 + (1j-1)*q*q/4 + b*b/2
    return lam, c, k0
def phi(z, b):
    lam, c, k0 = coeffs(b); return lam*z*z/2 - c*z + k0 - g*np.exp(-z*z/2)
def phi1(z, b):
    lam, c, k0 = coeffs(b); return lam*z - c + g*z*np.exp(-z*z/2)
def phi2(z, b):
    lam, c, k0 = coeffs(b); return lam + g*(1-z*z)*np.exp(-z*z/2)
def exact(b, hbar):
    lam, c, k0 = coeffs(b); mp.mp.dps = int(50 + 1.0*abs(g)/hbar)
    lam, c, k0, h = mp.mpc(lam), mp.mpc(c), mp.mpc(k0), mp.mpf(hbar); tot = mp.mpc(0)
    for p in range(int(60 + 8*abs(g)/hbar)):
        A = p - 1j*lam/h
        tot += (-1j*g/h)**p/mp.factorial(p)*mp.sqrt(2*mp.pi/A)*mp.e**(-(c/h)**2/(2*A))
    return complex(tot*mp.e**(1j*k0/h)), tot*mp.e**(1j*k0/h)
def saddles(b, box=4.0, ybox=3.0, step=0.25):
    found = []
    for re in np.arange(-box, box+1e-9, step):
        for im in np.arange(-ybox, ybox+1e-9, step):
            z = complex(re, im)
            try:
                for _ in range(60):
                    dz = phi1(z, b)/phi2(z, b); z -= dz
                    if abs(z) > 50: break
                if abs(phi1(z, b)) < 1e-10 and abs(z.real) <= box+1 and abs(z.imag) <= ybox+1 and all(abs(z-w) > 1e-6 for w in found): found.append(z)
            except (OverflowError, ZeroDivisionError): pass
    return sorted(found, key=lambda z: phi(z, b).imag)
def lead(z, b, hbar):
    p2 = phi2(z, b); th = np.pi/4 - np.angle(p2)/2
    if np.cos(th) < 0: th += np.pi
    return np.sqrt(2*np.pi*hbar/abs(p2))*np.exp(1j*th)*np.exp(1j*phi(z, b)/hbar)
if __name__ == "__main__":
    xk = a + p0; pk = p0 + g*xk*np.exp(-xk*xk/2)
    print(f"classical: packet centre at kick x_k = {xk}, momentum after kick {pk:.4f}, classical final position b_cl = {xk+pk:.4f}")
    for b in (-1.8, -0.5, 0.5, 1.0, 2.0):
        S = saddles(b)
        print(f"\nb = {b}: saddles with smallest Im phi (z, Im phi, |phi''|):")
        for z in S[:7]: print(f"    z = {z:.5f}  Im phi = {phi(z,b).imag:+.5f}  |phi''| = {abs(phi2(z,b)):.3f}")
        pos = [z for z in S if phi(z, b).imag > -1e-9][:5]
        for hbar in (0.1, 0.05, 0.025):
            I, _ = exact(b, hbar)
            best = None
            for k in range(1, min(4, len(pos))+1):
                for sub in itertools.combinations(range(len(pos)), k):
                    err = abs(I - sum(lead(pos[i], b, hbar) for i in sub))
                    if best is None or err < best[0]: best = (err, sub)
            single = abs(I - lead(pos[0], b, hbar))
            print(f"  hbar={hbar}: |I| = {abs(I):.4e}   error with lowest saddle only = {single:.3e}   best subset {best[1]} error = {best[0]:.3e}")
