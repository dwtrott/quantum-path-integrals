"""Interval certification of the selection certificate (Theorem S, correspondence 008) on Nova's polygonal contour
(correspondence 009).  All quantities are enclosed with Arb ball arithmetic through python-flint.
Parameters are exact rationals: lam = 1, g = -2, c = 553223/1000000, theta_ray = pi/8, kappa = 1/4,
local half-lengths rho = 51/100 (real saddle) and 21/100 (complex saddle),
external vertices -5/2 - 13/10 i, 3/5 - 13/10 i, 27/10 + 7/5 i.
Output: certified upper bounds on the certificate terms and a certified lower bound on the complex saddle term.
The code itself is unreviewed."""
from flint import arb, acb, fmpq, ctx
ctx.prec = 120
Q = lambda p, q=1: arb(fmpq(p, q))
C_ = Q(553223, 1000000); G_ = Q(-2); LAM = Q(1); KAPPA = Q(1, 4)
PI = arb.pi(); THETA = PI/8
one = acb(1)
def phi(z):  return LAM*z*z/2 - C_*z - G_*(-z*z/2).exp()
def phi1(z): return LAM*z - C_ + G_*z*(-z*z/2).exp()
def phi2(z): return LAM + G_*(one - z*z)*(-z*z/2).exp()
def phi3(z): return G_*(-z*z/2).exp()*(z*z*z - 3*z)
def phi4(z): return G_*(-z*z/2).exp()*(-(z*z*z*z) + 6*z*z - 3)
def up(x):  return arb(x).upper()      # arb enclosing an upper bound
def lo(x):  return arb(x).lower()
def ball(mid_re, mid_im, rad):  return acb(arb(mid_re, rad), arb(mid_im, rad))
def interval(a, b):  # arb ball covering [a, b], a <= b arbs (points)
    m = (a + b)/2; r = (b - a)/2
    return m + arb(0, float(up(r)) * (1 + 1e-9) + 1e-30)

def certify_root(approx_re, approx_im, rad):
    """interval Newton for the holomorphic phi1: unique zero in the box, returned as a refined enclosure."""
    Z = ball(approx_re, approx_im, rad)
    m = acb(arb(approx_re), arb(approx_im))
    d = phi2(Z)
    assert not d.contains(0), "derivative enclosure contains zero"
    N = m - phi1(m)/d
    inside = (lo(N.real) > lo(Z.real)) and (up(N.real) < up(Z.real)) and (lo(N.imag) > lo(Z.imag)) and (up(N.imag) < up(Z.imag))
    assert inside, "Newton image is not inside the box"
    return N          # encloses the unique zero in Z

# --- saddles (approximations from correspondence 009, boxes of radius 1e-9)
s_r = certify_root(1.514959245353315, 0.0, 1e-9)
s_c = certify_root(-0.634036499756894, -0.327421328952586, 1e-9)
print("real saddle enclosure   :", s_r.real.str(20), " imag radius", s_r.imag.str(5))
print("complex saddle enclosure:", s_c.real.str(20), s_c.imag.str(20))
# the real saddle is real: phi1 is real on R and changes sign on a real interval inside the box
xl, xh = Q(1514959245, 10**9), Q(1514959246, 10**9)
assert phi1(acb(xl)).real < 0 and phi1(acb(xh)).real > 0
s_r = acb(s_r.real, arb(0))

def local_data(s, rho, npanel=400):
    p2 = phi2(s); a = abs(p2); th = PI/4 - p2.arg()/2
    u = acb(th.cos(), th.sin())
    if not (th.cos() > 0): raise SystemExit("orientation not certified")
    M3 = arb(0); M4 = arb(0)
    for j in range(npanel):
        t = interval(-rho + 2*rho*j/npanel, -rho + 2*rho*(j+1)/npanel)
        z = s + u*t
        M3 = M3.max(up(abs(phi3(z)))); M4 = M4.max(up(abs(phi4(z))))
    assert rho*M3 < 3*KAPPA*a, "radius condition of Lemma L not certified"
    return dict(s=s, a=a, th=th, u=u, rho=rho, M3=M3, M4=M4, zm=s - u*rho, zp=s + u*rho, imphi=phi(s).imag)

Lr = local_data(s_r, Q(51, 100)); Lc = local_data(s_c, Q(21, 100))
for n, d in (("x_r", Lr), ("z_1", Lc)):
    print(f"{n}: a in {d['a'].str(12)}, theta in {d['th'].str(10)}, M3 <= {d['M3'].str(8)}, M4 <= {d['M4'].str(8)}, Im phi(s) in {d['imphi'].str(12)}")
assert Lr['imphi'].contains(0) or abs(Lr['imphi']) < 1e-30

V1 = acb(Q(-5, 2), Q(-13, 10)); V2 = acb(Q(3, 5), Q(-13, 10)); V3 = acb(Q(27, 10), Q(7, 5))
connectors = [("V1 -> z_1 segment", V1, Lc['zm']), ("z_1 segment -> V2", Lc['zp'], V2), ("V2 -> x_r segment", V2, Lr['zm']), ("x_r segment -> V3", Lr['zp'], V3)]

def connector_panels(A, B, tol=0.002, maxdepth=40):
    """adaptive panels; returns list of (length_upper, m_lower) with Im phi >= m_lower on the panel."""
    L = up(abs(B - A)); out = []
    stack = [(Q(0), Q(1), 0)]
    while stack:
        u0, u1, depth = stack.pop()
        t = interval(u0, u1); z = A*(1 - t) + B*t
        im = phi(z).imag
        width = float(up(im)) - float(lo(im))
        if width > tol and depth < maxdepth: 
            um = (u0 + u1)/2; stack += [(u0, um, depth+1), (um, u1, depth+1)]
        else:
            out.append((L*(u1 - u0), lo(im)))
    return out

panels = {name: connector_panels(A, B) for name, A, B in connectors}
mu = None
for name, ps in panels.items():
    m = min(float(p[1]) for p in ps); mu = m if mu is None else min(mu, m)
    print(f"connector {name}: {len(ps)} panels, certified min Im phi >= {m:.6f}")
print(f"certified mu (connectors) >= {mu:.6f}")

def ray_constants(E, sign):
    x0, y0 = E.real, E.imag; ct, st = THETA.cos(), THETA.sin()
    q0 = (x0 - C_)*y0; q1 = sign*((x0 - C_)*st + y0*ct); q2 = st*ct
    grow = sign*(x0*ct - y0*st)            # d/dt of (x^2-y^2)/2 at t=0 is grow; must be >= 0
    assert q1 > 0 and grow > 0, "ray monotonicity conditions not certified"
    G0 = ((y0*y0 - x0*x0)/2).exp()
    m0 = q0 - 2*G0                         # Im phi >= m0 + q2 t^2 on the ray (|g| = 2, lam = 1)
    return lo(m0), q2
rays = {"left": ray_constants(V1, -1), "right": ray_constants(V3, +1)}
for k, (m0, q2) in rays.items(): print(f"{k} ray: Im phi >= {float(m0):.4f} + {float(lo(q2)):.4f} t^2")

print("\n hbar      | certified bound B | local x_r | local z_1 | tail T | certified lower bound on |z_1 term| | B < |z_1 term| certified")
for hb in (Q(1, 20), Q(1, 40), Q(1, 80), Q(1, 160)):
    def local_err(d):
        a = d['a']; pref = (-d['imphi']/hb).exp()*(2*PI*hb/a).sqrt()
        r = (-a*d['rho']**2/(2*hb)).exp() + hb*(d['M4']/(8*a*a) + Q(5, 24)*d['M3']**2/((1 - KAPPA)**Q(7, 2)*a**3))
        return up(pref*r)
    Er, Ec = local_err(Lr), local_err(Lc)
    T = arb(0)
    for ps in panels.values():
        for length, m in ps: T += length*(-m/hb).exp()
    for m0, q2 in rays.values(): T += (-m0/hb).exp()*(PI*hb/q2).sqrt()/2
    T = up(T); B = up(Er + Ec + T)
    zterm = lo((-Lc['imphi']/hb).exp()*(2*PI*hb/Lc['a']).sqrt())
    print(f" {float(hb):<9} | {float(B):.6e} | {float(Er):.3e} | {float(Ec):.3e} | {float(T):.3e} | {float(zterm):.6e} | {bool(B < zterm)} | certified lower bound on the error of the real-saddle-only leading term: {float(lo(zterm - B)):.4e}")
