"""Checks for correspondence 003 (Nova's tree-cumulant bound). Floating point / mpmath evidence only.
Model: harmonic oscillator + eps*V0*exp(-x^2/(2 sigma^2)), trapezoid slicing, endpoint-stripped ratio Z.
Z = L(exp(-i eps F/hbar)), F = d*sum_j W(x_j), L = normalized Fresnel functional (mean x*, covariance i*hbar*G).
"""
import itertools, math
import mpmath as mp
mp.mp.dps = 60
m = w = V0 = mp.mpf(1)

def setup(N, T, a, b):
    d = mp.mpf(T)/N; n = N-1
    th = mp.acos(1-w*w*d*d/2)
    A = mp.matrix(n, n)
    for j in range(n):
        A[j, j] = (m/d)*(2-w*w*d*d)
        if j+1 < n: A[j, j+1] = A[j+1, j] = -(m/d)
    c = mp.matrix(n, 1); c[0] += (m/d)*a; c[n-1] += (m/d)*b
    G = A**-1; xs = G*c
    return d, n, th, A, c, G, xs

def gauss_moment(pvec, G, xs, hbar, sig):
    """L( prod_j W(x_j)^{p_j} ) / V0^{sum p}, branch continuous from hbar = 0."""
    n = len(pvec)
    idx = [j for j in range(n) if pvec[j] > 0]
    if not idx: return mp.mpf(1)
    k = len(idx)
    sq = [mp.sqrt(mp.mpf(pvec[j]))/sig for j in idx]
    Mm = mp.matrix(k, k)
    for a_ in range(k):
        for b_ in range(k):
            Mm[a_, b_] = sq[a_]*G[idx[a_], idx[b_]]*sq[b_]
    y = mp.matrix(k, 1)
    for a_ in range(k): y[a_] = sq[a_]*xs[idx[a_]]
    if k == 1:
        ev = [Mm[0, 0]]
    else:
        ev = mp.eigsy(Mm, eigvals_only=True)
    pref = mp.mpf(1)
    for g in ev: pref *= 1/mp.sqrt(1+1j*hbar*g)          # principal root, Re > 0
    K = mp.eye(k) + 1j*hbar*Mm
    v = mp.lu_solve(K, y)
    return pref*mp.e**(-(y.T*v)[0]/2)

def Z_coeffs(N, T, a, b, hbar, sig, P):
    """Taylor coefficients a_p of Z in eps, p = 0..P (exact Fourier/Dyson series, all slices)."""
    d, n, th, A, c, G, xs = setup(N, T, a, b)
    coef = [mp.mpc(0)]*(P+1)
    for p in range(P+1):
        tot = mp.mpc(0)
        for pv in itertools.product(range(p+1), repeat=n):
            if sum(pv) != p: continue
            wgt = mp.mpf(1)
            for q in pv: wgt /= mp.factorial(q)
            tot += wgt*gauss_moment(pv, G, xs, hbar, sig)
        coef[p] = (-1j*d*V0/hbar)**p*tot
    return coef, (d, n, th, A, c, G, xs)

def log_series(acoef):
    P = len(acoef)-1; c = [mp.mpc(0)]*(P+1)
    for p in range(1, P+1):
        c[p] = acoef[p] - sum(k*c[k]*acoef[p-k] for k in range(1, p))/p
    return c

def saddle(N, T, a, b, eps, hbar, sig):
    d, n, th, A, c, G, xs = setup(N, T, a, b)
    W = lambda x: V0*mp.e**(-x*x/(2*sig**2)); W1 = lambda x: -x/sig**2*W(x); W2 = lambda x: (x*x/sig**4-1/sig**2)*W(x)
    x = xs.copy()
    for _ in range(200):
        g = A*x - c - eps*d*mp.matrix([W1(x[j]) for j in range(n)])
        H = A.copy()
        for j in range(n): H[j, j] -= eps*d*W2(x[j])
        x = x - mp.lu_solve(H, g)
    S = lambda z, e: ((z.T*A*z)[0]/2 - (c.T*z)[0] - e*d*sum(W(z[j]) for j in range(n)))
    H = A.copy()
    for j in range(n): H[j, j] -= eps*d*W2(x[j])
    logsc = 1j*(S(x, eps)-S(xs, 0))/hbar - mp.log(abs(mp.det(H)/mp.det(A)))/2
    return logsc

def nova_bound(p, d_eff, M, hbar, sig):
    """|c_p| <= (r^2/(hbar M)) (d_eff L_r M/r^2)^p (p-2)!/p! C(3p-3,p-2), minimized over r; d_eff = n*d."""
    best = None
    for t in [mp.mpf(i)/20 for i in range(4, 80)]:
        r = sig*t; Lr = 2*V0*mp.e**(t*t/2)*mp.ncdf(t)
        val = (r*r/(hbar*M))*(d_eff*Lr*M/(r*r))**p*mp.factorial(p-2)/mp.factorial(p)*mp.binomial(3*p-3, p-2)
        best = val if best is None or val < best else best
    return best

print("== Pruefer identity: sum over labelled trees of prod deg! = (p-2)! C(3p-3,p-2)")
for p in range(2, 8):
    tot = 0
    for seq in itertools.product(range(p), repeat=p-2):
        deg = [1]*p
        for s in seq: deg[s] += 1
        tot += math.prod(math.factorial(x) for x in deg)
    print(p, tot, math.factorial(p-2)*math.comb(3*p-3, p-2))

print("== (1) log-series coefficients |c_p| vs tree bound; n*d used for the site sum; M = max|G_jl|")
for (N, T, lab) in ((2, 1.0, "N=2 before caustic"), (2, 3.2, "N=2 beyond discrete caustic (w d=1.6)"), (3, 1.0, "N=3"), (4, 4.0, "N=4 beyond first caustic")):
    for hbar in (1.0, 0.1):
        sig = mp.mpf(1); P = 9 if N <= 3 else 7
        ac, (d, n, th, A, c, G, xs) = Z_coeffs(N, T, 0.3, -0.7, hbar, sig, P)
        cl = log_series(ac); M = max(abs(G[i, j]) for i in range(n) for j in range(n))
        rat = [mp.nstr(abs(cl[p])/nova_bound(p, n*d, M, hbar, sig), 3) for p in range(2, P+1)]
        print(lab, "hbar=", hbar, "nu=", int(mp.floor(N*th/mp.pi)), " max_p |c_p|/bound_p:", max(float(x) for x in rat), " hbar*|c_2..4|:", [mp.nstr(hbar*abs(cl[p]), 4) for p in (2, 3, 4)])

print("== (2) radius of convergence of log Z (nearest zero of Z in complex eps), N=2, vs guaranteed radius")
for T in (1.0, 3.2):
    for hbar in (1.0, 0.3, 0.1):
        sig = mp.mpf(1); P = 46
        ac, (d, n, th, A, c, G, xs) = Z_coeffs(2, T, 0.3, -0.7, hbar, sig, P)
        cl = log_series(ac); M = abs(G[0, 0])
        est = [abs(cl[p])**(-mp.mpf(1)/p) for p in (40, 43, 46)]
        guar = max((sig*t)**2/(mp.mpf(27)/4*d*2*V0*mp.e**(t*t/2)*mp.ncdf(t)*M) for t in [mp.mpf(i)/20 for i in range(4, 80)])
        uniq = sig**2/(d*M*V0)
        print(f"T={T} hbar={hbar}: root-test radius ~ {[mp.nstr(e,4) for e in est]}  tree-guaranteed radius {mp.nstr(guar,4)}  Lemma-D (max norm) threshold {mp.nstr(uniq,4)}")

print("== (3) semiclassical test: log Z (series, branch from eps=0) vs saddle prediction; abs error and ratio error")
for (N, T, eps) in ((2, 1.0, 0.3), (3, 1.0, 0.3), (4, 4.0, 0.02)):
    for hbar in (1.0, 0.3, 0.1, 0.03):
        if N >= 3 and hbar < 0.09: continue
        sig = mp.mpf(1)
        tau = eps*(T/N)/hbar*(N-1)
        P = int(min(60, 14+4*tau)) if N == 2 else (int(min(22, 10+3*tau)) if N == 3 else 9)
        ac, _ = Z_coeffs(N, T, 0.3, -0.7, hbar, sig, P)
        Z = sum(ac[p]*mp.mpf(eps)**p for p in range(P+1))
        last = abs(ac[P]*mp.mpf(eps)**P)
        ls = saddle(N, T, 0.3, -0.7, mp.mpf(eps), hbar, sig)
        Zsc = mp.e**ls
        print(f"N={N} T={T} eps={eps} hbar={hbar}: |Z-Zsc|={mp.nstr(abs(Z-Zsc),4)}  |Z/Zsc-1|={mp.nstr(abs(Z/Zsc-1),4)}  ratio/hbar={mp.nstr(abs(Z/Zsc-1)/hbar,4)}  (series truncation ~{mp.nstr(last,2)})")
