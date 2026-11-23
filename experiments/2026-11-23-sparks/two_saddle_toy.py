"""Scoping experiment for the multi-saddle direction (N = 2 slice, one interior variable).
I(hbar) = lim_{eta->0} int exp(i*phi(x)/hbar) dx,  phi(x) = lam*x^2/2 - c*x - g*exp(-x^2/(2 sig^2)).
In the project's variables lam = (m/d)(2 - w^2 d^2), c = (m/d)(a+b), g = eps*d*V0.  g < 0 is an attractive Gaussian well.
Exact reference: Fourier/Dyson series (entire in g), evaluated with mpmath.  Evidence only, not certification.
"""
import mpmath as mp
mp.mp.dps = 120
lam = mp.mpf(1); sig = mp.mpf(1)

def phi(x, c, g):  return lam*x*x/2 - c*x - g*mp.e**(-x*x/(2*sig**2))
def phi1(x, c, g): return lam*x - c + g*x/sig**2*mp.e**(-x*x/(2*sig**2))
def phi2(x, c, g): return lam + g/sig**2*(1 - x*x/sig**2)*mp.e**(-x*x/(2*sig**2))

def exact(c, g, hbar, P=None):
    P = P or int(40 + 6*abs(g)/hbar)
    tot = mp.mpc(0)
    for p in range(P+1):
        a = p/sig**2 - 1j*lam/hbar
        tot += (-1j*g/hbar)**p/mp.factorial(p)*mp.sqrt(2*mp.pi/a)*mp.e**(-(c/hbar)**2/(2*a))
    return tot

def contrib(z, c, g, hbar):
    """Leading contribution of a nondegenerate saddle z, descent direction oriented left to right."""
    p2 = phi2(z, c, g)
    th = mp.pi/4 - mp.arg(p2)/2
    if mp.cos(th) < 0: th += mp.pi
    return mp.sqrt(2*mp.pi*hbar/abs(p2))*mp.e**(1j*th)*mp.e**(1j*phi(z, c, g)/hbar)

def real_saddles(c, g):
    xs = []; grid = [mp.mpf(i)/50 for i in range(-400, 401)]
    for u, v in zip(grid[:-1], grid[1:]):
        if phi1(u, c, g)*phi1(v, c, g) < 0:
            xs.append(mp.findroot(lambda x: phi1(x, c, g), (u, v), solver='anderson'))
    return xs

def complex_saddles(c, g, box=4.5, step=0.3):
    found = []
    re = -box
    while re <= box:
        im = -box
        while im <= box:
            try:
                z = mp.findroot(lambda x: phi1(x, c, g), mp.mpc(re, im), tol=1e-40, maxsteps=60)
                if abs(mp.im(z)) > 1e-12 and abs(z) < 6 and all(abs(z-w) > 1e-8 for w in found): found.append(z)
            except Exception: pass
            im += step
        re += step
    return found

g = mp.mpf(-2)
# fold: phi' = phi'' = 0 at x_f > 0
xf = mp.findroot(lambda x: phi2(x, 0, g), mp.mpf("0.6")); cf = lam*xf + g*xf/sig**2*mp.e**(-xf*xf/(2*sig**2))
print("g =", g, " cusp at c=0 for g < -lam*sig^2 =", -lam*sig**2, "; fold point x_f =", mp.nstr(xf, 8), " |c_fold| =", mp.nstr(abs(cf), 8))

print("\n== inside the cusp region: three real saddles")
for c in (mp.mpf('0.2'),):
    xs = real_saddles(c, g); print("c =", c, "real saddles:", [mp.nstr(x, 6) for x in xs], "phi'':", [mp.nstr(phi2(x, c, g), 4) for x in xs])
    for hbar in (mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('0.025')):
        I = exact(c, g, hbar); parts = [contrib(x, c, g, hbar) for x in xs]
        best1 = min(abs(I-p) for p in parts)
        print(f"  hbar={mp.nstr(hbar,3)}: |I|={mp.nstr(abs(I),5)}  best single-saddle error={mp.nstr(best1,4)}  three-saddle error={mp.nstr(abs(I-sum(parts)),4)}  /(hbar^1.5)={mp.nstr(abs(I-sum(parts))/hbar**1.5,4)}")

print("\n== just outside the fold: one real saddle plus complex saddles")
for dc in (mp.mpf('0.05'), mp.mpf('0.15')):
    c = abs(cf) + dc
    xs = real_saddles(c, g); zs = complex_saddles(c, g)
    zs = sorted(zs, key=lambda z: abs(mp.im(phi(z, c, g))))
    print("c = |c_fold| +", dc, " real saddles:", [mp.nstr(x, 6) for x in xs])
    print("  complex saddles found (|z|<6):", len(zs), " nearest five  z, Im phi:", [(mp.nstr(z, 5), mp.nstr(mp.im(phi(z, c, g)), 4)) for z in zs[:6]])
    zdec = [z for z in zs if mp.im(phi(z, c, g)) > 0]
    z1 = zdec[0]
    for hbar in (mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('0.025'), mp.mpf('0.0125')):
        I = exact(c, g, hbar); r = sum(contrib(x, c, g, hbar) for x in xs); cz = contrib(z1, c, g, hbar)
        print(f"  hbar={mp.nstr(hbar,3)}: |I|={mp.nstr(abs(I),5)}  real-only error={mp.nstr(abs(I-r),4)}  |complex term|={mp.nstr(abs(cz),4)}  real+complex error={mp.nstr(abs(I-r-cz),4)}  exp(-Im phi/hbar)={mp.nstr(mp.e**(-mp.im(phi(z1,c,g))/hbar),4)}")
