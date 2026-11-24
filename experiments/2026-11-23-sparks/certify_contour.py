"""General interval certification of the selection certificate (Theorem S) on a polygonal contour.
phi(z) = z^2/2 - c z - g exp(-z^2/2)  (lam = 1, sigma = 1), I = regulated integral of exp(i phi/hbar) over R.
Arb ball arithmetic through python-flint, 140 bits.  Every hypothesis of Lemma L and Theorem S is checked as an
inequality between balls; the script aborts if any check is not certified.  The code itself is unreviewed.

Usage: python certify_contour.py (three0 | three02 | two) [L | L2]
"""
import sys
from flint import arb, acb, fmpq, ctx
ctx.prec = 140
Q = lambda p, q=1: arb(fmpq(p, q))
PI = arb.pi(); THETA = PI/8; KAPPA = Q(1, 4); one = acb(1); LEMMA = "L"
up = lambda x: arb(x).upper(); lo = lambda x: arb(x).lower()

def make(c, g):
    phi  = lambda z: z*z/2 - c*z - g*(-z*z/2).exp()
    phi1 = lambda z: z - c + g*z*(-z*z/2).exp()
    phi2 = lambda z: one + g*(one - z*z)*(-z*z/2).exp()
    phi3 = lambda z: g*(-z*z/2).exp()*(z*z*z - 3*z)
    phi4 = lambda z: g*(-z*z/2).exp()*(-(z*z*z*z) + 6*z*z - 3)
    return phi, phi1, phi2, phi3, phi4

def interval(a, b):
    m = (a + b)/2; r = (b - a)/2
    return m + arb(0, float(up(r))*(1 + 1e-9) + 1e-40)

def certify_root(F, re, im, rad=1e-9, real=False):
    phi, phi1, phi2 = F[0], F[1], F[2]
    Z = acb(arb(re, rad), arb(im, rad)); m = acb(arb(re), arb(im))
    d = phi2(Z); assert not d.contains(0), "derivative enclosure contains zero"
    N = m - phi1(m)/d
    assert lo(N.real) > lo(Z.real) and up(N.real) < up(Z.real) and lo(N.imag) > lo(Z.imag) and up(N.imag) < up(Z.imag), "Newton image not inside box"
    if real:   # phi1 is real on R: a sign change on a real subinterval of the box puts the unique zero on R
        a_, b_ = arb(re) - arb(rad)/2, arb(re) + arb(rad)/2
        fa, fb = phi1(acb(a_)).real, phi1(acb(b_)).real
        assert (fa < 0 and fb > 0) or (fa > 0 and fb < 0), "no certified sign change"
        N = acb(N.real, arb(0))
    return N

def local_data(F, s, rho, npanel=600):
    phi, _, phi2, phi3, phi4 = F
    p2 = phi2(s); a = abs(p2); th = PI/4 - p2.arg()/2; u = acb(th.cos(), th.sin())
    assert th.cos() > 0, "orientation not certified"
    M3 = arb(0); M4 = arb(0)
    for j in range(npanel):
        t = interval(-rho + 2*rho*j/npanel, -rho + 2*rho*(j+1)/npanel); z = s + u*t
        M3 = M3.max(up(abs(phi3(z)))); M4 = M4.max(up(abs(phi4(z))))
    assert rho*M3 < 3*KAPPA*a, f"radius condition of Lemma L not certified (rho M3 = {float(up(rho*M3)):.4f}, 3 kappa a = {float(lo(3*KAPPA*a)):.4f})"
    return dict(s=s, a=a, th=th, u=u, rho=rho, M3=M3, M4=M4, zm=s - u*rho, zp=s + u*rho, imphi=phi(s).imag)

def connector(F, A, B, tol=0.002, maxdepth=44):
    phi = F[0]; L = up(abs(B - A)); out = []; stack = [(Q(0), Q(1), 0)]
    while stack:
        u0, u1, depth = stack.pop(); t = interval(u0, u1); im = phi(A*(1 - t) + B*t).imag
        if float(up(im)) - float(lo(im)) > tol and depth < maxdepth:
            um = (u0 + u1)/2; stack += [(u0, um, depth+1), (um, u1, depth+1)]
        else: out.append((L*(u1 - u0), lo(im)))
    return out

def ray(c, g, E, sign):
    x0, y0 = E.real, E.imag; ct, st = THETA.cos(), THETA.sin()
    q0 = (x0 - c)*y0; q1 = sign*((x0 - c)*st + y0*ct); q2 = st*ct; grow = sign*(x0*ct - y0*st)
    assert q1 > 0 and grow > 0, "ray monotonicity conditions not certified"
    m0 = q0 - abs(g)*((y0*y0 - x0*x0)/2).exp()
    return lo(m0), q2

def run(name, c, g, saddles, VL, VR, hbars, extra_vertices=None):
    """saddles: list of (re, im, real?, rho) ordered along the contour; extra_vertices: dict index -> vertex inserted after saddle index."""
    F = make(c, g); phi = F[0]
    print(f"=== {name}: c = {c.str(10)}, g = {g.str(4)}")
    loc = []
    for (re, im, real, rho) in saddles:
        s = acb(arb(0), arb(0)) if (re == 0 and im == 0) else certify_root(F, re, im, real=real)
        if re == 0 and im == 0: assert F[1](s).contains(0) and abs(F[1](s)) < 1e-60 and c == 0
        d = local_data(F, s, rho); loc.append(d)
        print(f" saddle {float(s.real.mid()):+.9f}{float(s.imag.mid()):+.9f}i: a = {float(d['a'].mid()):.6f}, theta = {float(d['th'].mid()):+.5f}, rho = {float(rho.mid()):.3f}, M3 <= {float(d['M3']):.4f}, M4 <= {float(d['M4']):.4f}, Im phi(s) in {d['imphi'].str(8)}")
    pts = [("VL", VL)]
    for i, d in enumerate(loc):
        pts.append(("in", d['zm'])); pts.append(("out", d['zp']))
        if extra_vertices and i in extra_vertices: pts.append(("V", extra_vertices[i]))
    pts.append(("VR", VR))
    panels = []; mu = None
    for (ka, A), (kb, B) in zip(pts[:-1], pts[1:]):
        if ka == "in" and kb == "out": continue        # saddle segment
        ps = connector(F, A, B); panels += ps
        mseg = min(float(p[1]) for p in ps); mu = mseg if mu is None else min(mu, mseg)
        print(f"  connector {ka}->{kb}: {len(ps)} panels, Im phi >= {mseg:.6f}")
    rays = [ray(c, g, VL, -1), ray(c, g, VR, +1)]
    print(f"  rays: Im phi >= {float(rays[0][0]):.4f} + {float(lo(rays[0][1])):.4f} t^2 (left), >= {float(rays[1][0]):.4f} + ... (right);  certified mu = {mu:.6f}")
    rows = []
    for hb in hbars:
        E = []; Lmod = []
        for d in loc:
            a = d['a']; pref = (-d['imphi']/hb).exp()*(2*PI*hb/a).sqrt()
            if LEMMA == "L":     # Lemma L of correspondence 008
                r = (-a*d['rho']**2/(2*hb)).exp() + hb*(d['M4']/(8*a*a) + Q(5, 24)*d['M3']**2/((1 - KAPPA)**Q(7, 2)*a**3))
            else:                # Lemma L2 of correspondence 015: sharp first-order constants plus an explicit hbar^(3/2) term
                r = ((-a*d['rho']**2/(2*hb)).exp() + hb*(d['M4']/(8*a*a) + Q(5, 24)*d['M3']**2/a**3)
                     + Q(32, 9)/(2*PI).sqrt()*d['M3']**3*hb**Q(3, 2)/((1 - KAPPA)**5*a**Q(9, 2)))
            E.append(up(pref*r)); Lmod.append(lo(pref))
        T = arb(0)
        for length, m in panels: T += length*(-m/hb).exp()
        for m0, q2 in rays: T += (-m0/hb).exp()*(PI*hb/q2).sqrt()/2
        T = up(T); B = up(sum(E, arb(0)) + T)
        rows.append((hb, B, E, T, Lmod))
        print(f"  hbar = 1/{int(round(1/float(hb)))}: certified B = {float(B):.6e}  (local {[f'{float(e):.3e}' for e in E]}, tail {float(T):.3e});  certified lower bounds on |L_s|: {[f'{float(l):.5e}' for l in Lmod]};  B < min|L_s|: {bool(B < min(Lmod, key=lambda v: float(v)))}")
    return loc, rows

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "three0"
    globals()['LEMMA'] = sys.argv[2] if len(sys.argv) > 2 else "L"
    hb = [Q(1, 20), Q(1, 40), Q(1, 80), Q(1, 160), Q(1, 320)]
    if which == "three0":      # symmetric well, three real saddles: 0 and +-sqrt(2 log 2)
        run("three real saddles, symmetric", Q(0), Q(-2),
            [(-1.1774100225154747, 0.0, True, Q(33, 100)), (0.0, 0.0, True, Q(33, 100)), (1.1774100225154747, 0.0, True, Q(33, 100))],
            acb(Q(-28, 10), Q(-1)), acb(Q(28, 10), Q(1)), hb)
    elif which == "three02":   # asymmetric, inside the cusp region
        import mpmath as mp      # approximations only; every root is certified below by interval Newton
        mp.mp.dps = 30
        f = lambda x: x - mp.mpf(1)/5 - 2*x*mp.e**(-x*x/2)
        r = [float(mp.findroot(f, x0)) for x0 in (-1.01, -0.21, 1.31)]
        run("three real saddles, c = 1/5", Q(1, 5), Q(-2),
            [(r[0], 0.0, True, Q(25, 100)), (r[1], 0.0, True, Q(25, 100)), (r[2], 0.0, True, Q(40, 100))],
            acb(Q(-28, 10), Q(-1)), acb(Q(28, 10), Q(1)), hb)
    else:                      # the two-saddle instance of correspondence 010, for regression
        run("one real and one complex saddle (correspondence 010)", Q(553223, 1000000), Q(-2),
            [(-0.634036499756894, -0.327421328952586, False, Q(21, 100)), (1.514959245353315, 0.0, True, Q(51, 100))],
            acb(Q(-5, 2), Q(-13, 10)), acb(Q(27, 10), Q(7, 5)), hb[:4], extra_vertices={0: acb(Q(3, 5), Q(-13, 10))})
