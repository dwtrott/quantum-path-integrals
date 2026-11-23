"""One-variable selection certificate (correspondence 008).  Evidence in floating point; NOT yet validated numerics.
phi(z) = lam z^2/2 - c z - g exp(-z^2/2)   (sigma = 1),   I = int_R exp(i phi/hbar) dx (regulated).
Builds an explicit contour through the real saddle x_r and the complex saddle z_1 (Im phi > 0), checks that the
contour integral reproduces the exact series, and evaluates the certificate bound of Theorem S."""
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
mp.mp.dps = 60
import sys
lam, g = 1.0, -2.0
kappa = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5   # parameter of the local lemma, 0 < kappa < 1
phi  = lambda z, c: lam*z*z/2 - c*z - g*np.exp(-z*z/2)
phi1 = lambda z, c: lam*z - c + g*z*np.exp(-z*z/2)
phi2 = lambda z, c: lam + g*(1 - z*z)*np.exp(-z*z/2)
phi3 = lambda z, c: g*np.exp(-z*z/2)*(z**3 - 3*z)
phi4 = lambda z, c: g*np.exp(-z*z/2)*(-z**4 + 6*z*z - 3)

def exact(c, hbar):
    mp.mp.dps = int(60 + 1.0*abs(g)/hbar)  # the series has cancellations of size exp(|g|/hbar)
    c, hbar = mp.mpf(c), mp.mpf(hbar); tot = mp.mpc(0)
    for p in range(int(60 + 8*abs(g)/hbar)):
        a = p - 1j*lam/hbar
        tot += (-1j*mp.mpf(g)/hbar)**p/mp.factorial(p)*mp.sqrt(2*mp.pi/a)*mp.e**(-(c/hbar)**2/(2*a))
    return complex(tot)

def newton(z, c):
    for _ in range(100): z = z - phi1(z, c)/phi2(z, c)
    return z

def descent(z0, c, Himag, smax=40.0):
    """steepest descent of Re(i phi) (Im phi increases), unit speed, until Im phi >= Himag."""
    def rhs(s, y):
        z = y[0] + 1j*y[1]; v = -np.conj(1j*phi1(z, c)); v /= abs(v); return [v.real, v.imag]
    ev = lambda s, y: phi(y[0] + 1j*y[1], c).imag - Himag
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(rhs, [0, smax], [z0.real, z0.imag], events=ev, rtol=1e-10, atol=1e-12, max_step=0.01)
    return sol.y[0] + 1j*sol.y[1]

def seg_int(f, pts):
    x, w = np.polynomial.legendre.leggauss(8); tot = 0
    for p, q in zip(pts[:-1], pts[1:]):
        zz = (p+q)/2 + (q-p)/2*x; tot += np.sum(w*f(zz))*(q-p)/2
    return tot

xf = brentq(lambda x: phi2(x, 0.0).real, 0.3, 1.0); cf = abs(lam*xf + g*xf*np.exp(-xf*xf/2))
c = cf + 0.15
xr = brentq(lambda x: phi1(x, c).real, 1.0, 3.0)
z1 = newton(-0.634 - 0.327j, c)
print(f"c = {c:.6f}; real saddle x_r = {xr:.6f}, phi''={phi2(xr,c).real:.4f}; complex saddle z_1 = {z1:.6f}, Im phi(z_1) = {phi(z1,c).imag:.6f}, phi''(z_1) = {phi2(z1,c):.4f}")

H = 4.0
sad = {}
for name, s in (("x_r", xr + 0j), ("z_1", z1)):
    p2 = phi2(s, c); a = abs(p2); th = np.pi/4 - np.angle(p2)/2
    if np.cos(th) < 0: th += np.pi
    rho = 0.0
    for cand in np.linspace(0.01, 0.8, 80):
        zz = s + np.linspace(-cand, cand, 801)*np.exp(1j*th)
        if cand <= 3*kappa*a/np.max(abs(phi3(zz, c))): rho = cand
        else: break
    t = np.linspace(-rho, rho, 801); zz = s + t*np.exp(1j*th)
    M3 = np.max(abs(phi3(zz, c))); M4 = np.max(abs(phi4(zz, c)))
    assert 0 < rho <= 3*kappa*a/M3, "segment too long for the local lemma"
    sad[name] = dict(s=s, a=a, th=th, rho=rho, M3=M3, M4=M4, zm=s - rho*np.exp(1j*th), zp=s + rho*np.exp(1j*th))
    print(f"{name}: a={a:.4f} theta={th:.4f} rho={rho:.4f} M3={M3:.4f} M4={M4:.4f}  Im phi at segment ends: {phi(sad[name]['zm'],c).imag:.4f}, {phi(sad[name]['zp'],c).imag:.4f}")

# tails: descent paths from the four segment ends
paths = {k: descent(sad[n][e], c, H) for k, (n, e) in dict(xr_plus=("x_r", "zp"), xr_minus=("x_r", "zm"), z1_plus=("z_1", "zp"), z1_minus=("z_1", "zm")).items()}
for k, p in paths.items(): print(f"  {k}: ends at {p[-1]:.4f}  (arg {np.degrees(np.angle(p[-1])):.1f} deg), Im phi = {phi(p[-1],c).imag:.3f}, length {np.sum(abs(np.diff(p))):.3f}")

# full contour
def ray(E, direction, tmax=7.0, n=1400): return E + np.linspace(0, tmax, n)*direction
left_ray  = ray(paths["z1_minus"][-1], -np.exp(1j*np.pi/8))[::-1]
right_ray = ray(paths["xr_plus"][-1],   np.exp(1j*np.pi/8))
conn = np.linspace(paths["z1_plus"][-1], paths["xr_minus"][-1], 200)
seg = lambda n: np.linspace(sad[n]["zm"], sad[n]["zp"], 400)
pieces = [("left ray", left_ray), ("z1_minus", paths["z1_minus"][::-1]), ("SEG z_1", seg("z_1")), ("z1_plus", paths["z1_plus"]),
          ("connector", conn), ("xr_minus", paths["xr_minus"][::-1]), ("SEG x_r", seg("x_r")), ("xr_plus", paths["xr_plus"]), ("right ray", right_ray)]
for name, p in pieces:
    if not name.startswith("SEG"): print(f"  min Im phi on {name}: {phi(p, c).imag.min():.4f}")
mu = min(phi(p, c).imag.min() for name, p in pieces if not name.startswith("SEG"))
print(f"mu = min Im phi off the saddle segments = {mu:.4f};  Im phi at saddles: x_r 0, z_1 {phi(z1,c).imag:.4f}")
for name, mono in (("left ray", left_ray[::-1]), ("right ray", right_ray)):
    print(f"  Im phi monotone increasing along {name}: {bool(np.all(np.diff(phi(mono, c).imag) > 0))}")

print("\n hbar | |I - contour integral| | actual error of two-saddle formula | certificate bound | local x_r | local z_1 | tail T | bound/actual")
for hbar in (0.05, 0.025, 0.0125, 0.00625):
    f = lambda z: np.exp(1j*phi(z, c)/hbar)
    Ic = sum(seg_int(f, p) for _, p in pieces); Iex = exact(c, hbar)
    L = {n: np.exp(1j*phi(d["s"], c)/hbar)*np.sqrt(2*np.pi*hbar/d["a"])*np.exp(1j*d["th"]) for n, d in sad.items()}
    E = {n: np.exp(-phi(d["s"], c).imag/hbar)*np.sqrt(2*np.pi*hbar/d["a"])*(np.exp(-d["a"]*d["rho"]**2/(2*hbar)) + hbar*(d["M4"]/(8*d["a"]**2) + 5/24*(1-kappa)**-3.5*d["M3"]**2/d["a"]**3)) for n, d in sad.items()}
    T = sum(np.sum(np.exp(-phi((p[:-1]+p[1:])/2, c).imag/hbar)*abs(np.diff(p))) for name, p in pieces if not name.startswith("SEG"))
    act = abs(Iex - sum(L.values())); B = sum(E.values()) + T
    # check of the local lemma itself on each segment
    chk = {n: (abs(seg_int(f, seg(n)) - L[n]), E[n]) for n in sad}
    print(f" {hbar:<7}| {abs(Iex-Ic):.2e} | {act:.3e} | {B:.3e} | {E['x_r']:.3e} | {E['z_1']:.3e} | {T:.3e} | {B/act:.1f}   local lemma check (actual, bound): x_r ({chk['x_r'][0]:.2e}, {chk['x_r'][1]:.2e}) z_1 ({chk['z_1'][0]:.2e}, {chk['z_1'][1]:.2e})")
