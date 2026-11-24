"""Interval certification of a two-variable selection certificate: two kicks, wavepacket in and out.
Claim certified for each hbar in the list:   | I - L0 | <= B_rel * |L0|,
where I = int_{R^2} exp(i Phi/hbar) dx, L0 = 2 pi hbar det(R) exp(i Phi(z_s)/hbar) is the leading term of one complex
saddle z_s, and R = (-i Phi''(z_s))^{-1/2}.  Arb ball arithmetic (python-flint).  The code is unreviewed.

Surface: z(u) = x_s + Re(R) u + i (y_s + clip_B(Im(R) u)), u in R^2.  Local square |u|_inf <= rho (exact tangent plane).
"""
import sys, time
import numpy as np
from flint import arb, acb, fmpq, ctx
ctx.prec = 100
Q = lambda p, q=1: arb(fmpq(p, q))
up = lambda x: arb(x).upper(); lo = lambda x: arb(x).lower()
I_ = acb(0, 1); PI = arb.pi()
A_, P0, G_ = Q(-2), Q(1), Q(3)                       # a, p0, g
B_T, PB = Q(2), Q(1)                                 # target wavepacket: position b, momentum pb
KAPPA = Q(1, 4); BCLIP = Q(6, 10); RHO = Q(45, 100)
c14 = (I_ - 1)/4
def Phi(z1, z2):
    return ((z1 - A_)*(z1 - A_)/2 + c14*(P0 - (z1 - A_))*(P0 - (z1 - A_)) + (z2 - z1)*(z2 - z1)/2
            + (z2 - B_T)*(z2 - B_T)/2 + c14*(PB + (z2 - B_T))*(PB + (z2 - B_T)) - G_*((-z1*z1/2).exp() + (-z2*z2/2).exp()))
def grad(z1, z2):
    return ((z1 - A_) - 2*c14*(P0 - (z1 - A_)) - (z2 - z1) + G_*z1*(-z1*z1/2).exp(),
            (z2 - z1) + (z2 - B_T) + 2*c14*(PB + (z2 - B_T)) + G_*z2*(-z2*z2/2).exp())
def hess(z1, z2):
    return ((2 + 2*c14 + G_*(1 - z1*z1)*(-z1*z1/2).exp(), acb(-1)), (acb(-1), 2 + 2*c14 + G_*(1 - z2*z2)*(-z2*z2/2).exp()))
def inv2(M):
    (a, b), (c, d) = M; det = a*d - b*c
    return ((d/det, -b/det), (-c/det, a/det)), det
def mul2(M, N): return tuple(tuple(sum((M[i][k]*N[k][j] for k in range(2)), acb(0)) for j in range(2)) for i in range(2))
def interval(a, b):
    m = (a + b)/2; r = (b - a)/2
    return m + arb(0, float(up(r))*(1 + 1e-9) + 1e-40)

# ---- 1. saddle: Krawczyk test in C^2
def krawczyk(z0, rad=1e-8):
    m = (acb(z0[0].real, z0[0].imag), acb(z0[1].real, z0[1].imag))
    Z = tuple(acb(arb(float(z.real), rad), arb(float(z.imag), rad)) for z in z0)
    C, _ = inv2(hess(*m)); f = grad(*m); H = hess(*Z)
    CH = mul2(C, H); E = ((1 - CH[0][0], -CH[0][1]), (-CH[1][0], 1 - CH[1][1]))
    dZ = (Z[0] - m[0], Z[1] - m[1])
    K = tuple(m[i] - (C[i][0]*f[0] + C[i][1]*f[1]) + (E[i][0]*dZ[0] + E[i][1]*dZ[1]) for i in range(2))
    for i in range(2):
        assert lo(K[i].real) > lo(Z[i].real) and up(K[i].real) < up(Z[i].real) and lo(K[i].imag) > lo(Z[i].imag) and up(K[i].imag) < up(Z[i].imag), "Krawczyk test failed"
    return K
z0 = np.array([-0.28412 - 0.06944j, 0.28412 + 0.06944j])
for _ in range(30):      # float Newton, approximation only
    x1, x2 = z0
    gf = np.array([(x1+2) - (1j-1)*(1-(x1+2))/2 - (x2-x1) + 3*x1*np.exp(-x1*x1/2), (x2-x1) + (x2-2) + (1j-1)*(1+(x2-2))/2 + 3*x2*np.exp(-x2*x2/2)])
    Hf = np.array([[1.5+0.5j + 1 - 1 + 1 + 3*(1-x1*x1)*np.exp(-x1*x1/2) - 1, -1], [-1, 1.5+0.5j + 3*(1-x2*x2)*np.exp(-x2*x2/2)]])
    Hf[0, 0] = 2 + (1j-1)/2 + 3*(1-x1*x1)*np.exp(-x1*x1/2); Hf[1, 1] = 2 + (1j-1)/2 + 3*(1-x2*x2)*np.exp(-x2*x2/2)
    z0 = z0 - np.linalg.solve(Hf, gf)
zs = krawczyk(z0)
Phis = Phi(*zs); ImPhis = Phis.imag
print(f"saddle certified: z_s = ({zs[0].real.str(12)} + {zs[0].imag.str(12)} i, {zs[1].real.str(12)} + {zs[1].imag.str(12)} i);  Im Phi_s in {ImPhis.str(12)}")

# ---- 2. R = (-i Phi''(z_s))^{-1/2} by the 2x2 closed form (either sign of the determinant root gives a valid R)
H = hess(*zs); M = tuple(tuple(-I_*H[i][j] for j in range(2)) for i in range(2))
detM = M[0][0]*M[1][1] - M[0][1]*M[1][0]; s = detM.sqrt(); t = (M[0][0] + M[1][1] + 2*s).sqrt()
sqrtM = (((M[0][0] + s)/t, M[0][1]/t), (M[1][0]/t, (M[1][1] + s)/t))
R, _ = inv2(sqrtM); detR = R[0][0]*R[1][1] - R[0][1]*R[1][0]
chk = mul2(mul2(R, M), R)
assert all(abs(chk[i][j] - (1 if i == j else 0)) < 1e-9 for i in range(2) for j in range(2)), "R^T (-i Phi'') R = I not verified"
ReR = tuple(tuple(R[i][j].real for j in range(2)) for i in range(2)); ImR = tuple(tuple(R[i][j].imag for j in range(2)) for i in range(2))
detReR = ReR[0][0]*ReR[1][1] - ReR[0][1]*ReR[1][0]
assert detReR > 0, "orientation (det Re R > 0) not certified"
rown = [ (abs(R[j][0])**2 + abs(R[j][1])**2).sqrt() for j in range(2)]
print(f"|det R| in {abs(detR).str(10)}, det Re R in {detReR.str(8)}, row norms {[x.str(6) for x in rown]}")

# ---- 3. local data on the square |u|_inf <= rho
def W3(z): return (-z*z/2).exp()*(3*z - z*z*z)
def W4(z): return (-z*z/2).exp()*(z*z*z*z - 6*z*z + 3)
n_loc = 24; s3 = [arb(0), arb(0)]; s4 = [arb(0), arb(0)]
for i in range(n_loc):
    for k in range(n_loc):
        u = (interval(-RHO + 2*RHO*i/n_loc, -RHO + 2*RHO*(i+1)/n_loc), interval(-RHO + 2*RHO*k/n_loc, -RHO + 2*RHO*(k+1)/n_loc))
        for j in range(2):
            z = zs[j] + R[j][0]*u[0] + R[j][1]*u[1]
            s3[j] = s3[j].max(up(abs(W3(z)))); s4[j] = s4[j].max(up(abs(W4(z))))
M3 = up(G_*(s3[0]*rown[0]**3 + s3[1]*rown[1]**3)); M4 = up(G_*(s4[0]*rown[0]**4 + s4[1]*rown[1]**4))
assert Q(2).sqrt()*RHO*M3 < 3*KAPPA, f"radius condition not certified: sqrt2 rho M3 = {float(up(Q(2).sqrt()*RHO*M3)):.4f}"
for j in range(2): assert RHO*(abs(ImR[j][0]) + abs(ImR[j][1])) < BCLIP, "local square not inside the unclipped region"
print(f"rho = {float(RHO)}, M3 <= {float(M3):.4f}, M4 <= {float(M4):.4f}, sqrt2*rho*M3 = {float(up(Q(2).sqrt()*RHO*M3)):.4f} < 3 kappa = 0.75")

# ---- 4. tail over the surface outside the local square
xs = (zs[0].real, zs[1].real); ys = (zs[0].imag, zs[1].imag)
def panel(u1, u2):
    """returns (lower bound of Im Phi - Im Phi_s, width, upper bound of |det dz/du|)"""
    x = (xs[0] + ReR[0][0]*u1 + ReR[0][1]*u2, xs[1] + ReR[1][0]*u1 + ReR[1][1]*u2)
    v = (ImR[0][0]*u1 + ImR[0][1]*u2, ImR[1][0]*u1 + ImR[1][1]*u2)
    y = []; D = []
    for j in range(2):
        if abs(v[j]) < BCLIP: y.append(ys[j] + v[j]); D.append(arb(1))
        elif v[j] > BCLIP: y.append(ys[j] + BCLIP); D.append(arb(0))
        elif v[j] < -BCLIP: y.append(ys[j] - BCLIP); D.append(arb(0))
        else: y.append(ys[j] + v[j].min(BCLIP).max(-BCLIP)); D.append(arb(Q(1, 2), 0.5))
    z = (acb(x[0], y[0]), acb(x[1], y[1]))
    im = Phi(*z).imag - ImPhis
    J = (acb(ReR[0][0], D[0]*ImR[0][0])*acb(ReR[1][1], D[1]*ImR[1][1]) - acb(ReR[0][1], D[0]*ImR[0][1])*acb(ReR[1][0], D[1]*ImR[1][0]))
    return float(lo(im)), float(up(im)) - float(lo(im)), float(up(abs(J)))
HB = [Q(1, 20), Q(1, 40), Q(1, 80), Q(1, 160)]; hmin = float(HB[-1])
UMAX = 48.0; CUT = 1.3
leaves = []           # (area, dm_lower, Jmax)
stack = [(-UMAX, UMAX, -UMAX, UMAX)]; t0 = time.time(); rho_f = float(RHO)
while stack:
    a1, b1, a2, b2 = stack.pop()
    if a1 >= -rho_f and b1 <= rho_f and a2 >= -rho_f and b2 <= rho_f: continue     # inside the local square
    dm, wdt, J = panel(interval(Q(a1), Q(b1)) if False else arb((a1+b1)/2, (b1-a1)/2*(1+1e-12)), arb((a2+b2)/2, (b2-a2)/2*(1+1e-12)))
    straddle = not (b1 <= -rho_f or a1 >= rho_f or b2 <= -rho_f or a2 >= rho_f)   # overlaps the local square
    tol = 0.025*max(dm, 12*hmin) if dm > 0 else 0.0
    small = (b1 - a1) < 2e-4
    if (not straddle) and (dm > CUT or wdt <= tol or small):
        leaves.append(((b1-a1)*(b2-a2), dm, J))
    elif straddle and small:
        leaves.append(((b1-a1)*(b2-a2), dm, J))          # tiny boundary panel, counted fully
    else:
        m1, m2 = (a1+b1)/2, (a2+b2)/2
        stack += [(a1, m1, a2, m2), (m1, b1, a2, m2), (a1, m1, m2, b2), (m1, b1, m2, b2)]
print(f"{len(leaves)} panels in {time.time()-t0:.0f} s; smallest certified (Im Phi - Im Phi_s) on the surface outside the square: {min(l[1] for l in leaves):.5f}")
assert min(l[1] for l in leaves) > 0

# ---- 5. far region |u|_inf > UMAX: analytic Gaussian bound in x
Ymax = up((abs(ys[0]).max(abs(ys[1]))) + BCLIP); Lc = Q(5, 2)*Ymax
q1, q2 = A_ + P0, B_T - PB
CQ = (abs(P0) + abs(PB) + 2*abs(q2 - q1))*Ymax + Ymax*Ymax/2 + 2*G_*(Ymax*Ymax/2).exp()
# the u-box contains the x-box |x - q|_inf <= X when UMAX >= ||Re(R)^{-1}||_inf (X + |x_s - q|_inf)
ReRinv, _ = inv2(tuple(tuple(acb(ReR[i][j]) for j in range(2)) for i in range(2)))
ninf = max(float(up(abs(ReRinv[i][0]) + abs(ReRinv[i][1]))) for i in range(2))
off = max(float(up(abs(xs[0] - q1))), float(up(abs(xs[1] - q2))))
X = Q(int((UMAX/ninf - off)*1000), 1000); assert X > 2*Lc + 1
Jx = up((1 + abs(ImR[0][0])/1)*0 + 1)   # placeholder, replaced below
# |det(I + i D K)| with K = Im R (Re R)^{-1}, D in [0,1]^2: bounded by (1 + ||K||_F)^2
Kf = [[ImR[i][0]*ReRinv[0][j].real + ImR[i][1]*ReRinv[1][j].real for j in range(2)] for i in range(2)]
Jx = up((1 + (Kf[0][0]**2 + Kf[0][1]**2 + Kf[1][0]**2 + Kf[1][1]**2).sqrt())**2)

print(f"far region: |x - q|_inf > {float(X):.3f}, Im Phi >= sum s_j^2/4 - {float(Lc):.3f} sum|s_j| - {float(up(CQ)):.3f}")
print("\n hbar   | certified relative bound on |I - L0|/|L0| | local part | tail part | far part | |L0| in")
for hb in HB:
    r_loc = 2*(RHO/(2*hb).sqrt()).erfc() + hb*(M4/3 + Q(2, 3)*M3*M3/((1 - KAPPA)**4))
    L0abs = 2*PI*hb*abs(detR)*(-ImPhis/hb).exp()
    T = arb(0); hbf = float(hb)
    for area, dm, J in leaves: T += arb(area)*arb(J)*(-arb(dm)/hb).exp()
    T_rel = T/(2*PI*hb*abs(detR))                      # relative to |L0| (the factor exp(-Im Phi_s/hbar) is already removed)
    Gfull = 2*(Lc*Lc/hb).exp()*(4*PI*hb).sqrt(); Gtail = 2*(Lc*Lc/hb).exp()*(2*hb/(X - 2*Lc))*(-(X - 2*Lc)**2/(4*hb)).exp()
    far = (CQ/hb).exp()*2*Gtail*Gfull*Jx
    far_rel = far/L0abs
    tot = up(r_loc + T_rel + far_rel)
    print(f" 1/{int(round(1/hbf)):<4} | {float(tot):.5f} | {float(up(r_loc)):.5f} | {float(up(T_rel)):.3e} | {float(up(far_rel)):.1e} | {L0abs.str(6)}")
