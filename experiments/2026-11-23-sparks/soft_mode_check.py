"""Independent test of Nova's frozen-stiff soft-mode reduction (correspondence 011) for N = 3, 4.
J_N is computed from the exact Fourier series of the full (N-1)-dimensional integral (closed-form complex Gaussian
moments), J_soft by the same series with the stiff covariance set to zero, and the first covariance correction by a
derivative in a scaling parameter of the stiff covariance.  mpmath, 50 digits.  Evidence, not certification."""
import itertools, mpmath as mp
mp.mp.dps = 50
def series(Gmat, n, tau, hbar, P):
    tot = mp.mpc(0)
    for p in range(P+1):
        for pv in itertools.product(range(p+1), repeat=n):
            if sum(pv) != p: continue
            idx = [j for j in range(n) if pv[j] > 0]
            wgt = mp.mpf(1)
            for q in pv: wgt /= mp.factorial(q)
            if not idx: tot += wgt*(1j*tau)**p; continue
            k = len(idx); sq = [mp.sqrt(mp.mpf(pv[j])) for j in idx]
            Mm = mp.matrix(k, k)
            for a_ in range(k):
                for b_ in range(k): Mm[a_, b_] = sq[a_]*Gmat[idx[a_], idx[b_]]*sq[b_]
            det = mp.det(mp.eye(k) + 1j*hbar*Mm)
            # branch continuous from hbar = 0: product of principal roots over eigenvalues of the real symmetric Mm
            ev = [Mm[0, 0]] if k == 1 else mp.eigsy(Mm, eigvals_only=True)
            pref = mp.mpf(1)
            for g in ev: pref *= 1/mp.sqrt(1 + 1j*hbar*g)
            tot += wgt*(1j*tau)**p*pref
    return tot
def run(N, hbar):
    U = hbar; n = N-1
    T = 2*N*mp.sin(mp.pi/(2*N))/mp.sqrt(1+U/2); d = T/N; tau = U*d/hbar
    A = mp.matrix(n, n)
    for j in range(n):
        A[j, j] = (2 - d*d)/d
        if j+1 < n: A[j, j+1] = A[j+1, j] = -1/d
    v = mp.matrix([mp.sin((j+1)*mp.pi/N) for j in range(n)]); vv = (v.T*v)[0]
    lam1 = (v.T*A*v)[0]/vv; a = lam1*vv
    G = A**-1; Gsoft = (v*v.T)/(lam1*vv); Gperp = G - Gsoft
    M = max(abs(Gperp[i, j]) for i in range(n) for j in range(n))
    P = 34 if n == 2 else 26
    norm = mp.e**(1j*tau)*mp.sqrt(2*mp.pi*hbar/a)*mp.e**(1j*mp.pi/4)
    J = norm*series(G, n, tau, hbar, P)
    Js = norm*series(Gsoft, n, tau, hbar, P)
    eps = mp.mpf('1e-12')
    dJ = norm*(series(Gsoft + eps*Gperp, n, tau, hbar, P) - series(Gsoft - eps*Gperp, n, tau, hbar, P))/(2*eps)
    J1 = Js + dJ
    # Nova's bound with B = 1/4
    B = mp.mpf(1)/4; L0 = 2*mp.e**(B*B/2)*mp.ncdf(B); L1 = B*L0 + mp.sqrt(2/mp.pi); L2 = (1+B*B)*L0 + B*mp.sqrt(2/mp.pi)
    t = n*tau; D = hbar*M/2*mp.e**(t*L0)*(t*L2 + t*t*L1*L1); bound = (2*B + 2*hbar/(a*B))*D
    return dict(T=T, lam1_over_d=lam1/d, a_over_h=a/hbar, M=M, J=J, Js=Js, err0=abs(J-Js), err1=abs(J-J1), bound=bound)
print(" N  hbar=U   | J_N (exact series)            | |J-J_soft|  | Nova bound  | 6000 h | |J-J_1|   | M      (3T/8)")
for N in (3, 4):
    for hbar in (mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4')):
        r = run(N, hbar)
        print(f" {N}  {mp.nstr(hbar,3):8} | {mp.nstr(r['J'],9):29} | {mp.nstr(r['err0'],4):10} | {mp.nstr(r['bound'],4):11} | {mp.nstr(6000*hbar,3):6} | {mp.nstr(r['err1'],4):8} | {mp.nstr(r['M'],4)} ({mp.nstr(3*r['T']/8,4)})   check lam1/d={mp.nstr(r['lam1_over_d'],6)}, a/h={mp.nstr(r['a_over_h'],6)}, T/4={mp.nstr(r['T']/4,6)}")
print("Nova's table value for N=3 (U=h=1e-7 not stated; J_soft at h->0):  -4.619346312 + 0.087075125 i")
r = run(3, mp.mpf('1e-7')); print(" N=3, h=U=1e-7: J_soft =", mp.nstr(r['Js'], 10), " J =", mp.nstr(r['J'], 10), " |J-J_soft| =", mp.nstr(r['err0'], 4))
