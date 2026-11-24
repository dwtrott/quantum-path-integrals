"""Nova: spectral truncation benchmark, zero endpoints, positive harmonic Hessian.
Run: OPENBLAS_NUM_THREADS=1 python pa_benchmark.py
Floating point only. Analytic series tail and reduction bound exclude rounding.
"""
import itertools, math, platform, time
import numpy as np
import scipy
from scipy.special import roots_hermite, gammaln


def covariance(N, T):
    k = np.arange(1, N)
    d = T/N
    lam = (4*np.sin(np.pi*k/(2*N))**2-d*d)/d
    V = np.sqrt(2/N)*np.sin(np.pi*np.outer(k,k)/N)
    assert np.all(lam > 0)
    return V, lam, (V/lam)@V.T


def compositions(p, n):
    # stars and bars
    for b in itertools.combinations(range(p+n-1), n-1):
        yield np.diff((-1, *b, p+n-1))-1


def reference(G, tau, h, P=14):
    n = len(G)
    ans = 0j
    count = 0
    for p in range(P+1):
        a = np.asarray(list(compositions(p,n)))
        s = np.sqrt(a)
        ev = np.linalg.eigvalsh(s[:,:,None]*G[None,:,:]*s[:,None,:])
        moment = np.exp(-.5*np.log1p(1j*h*ev).sum(axis=1))
        ans += (1j*tau)**p*np.sum(moment*np.exp(-gammaln(a+1).sum(axis=1)))
        count += len(a)
    t = n*abs(tau)
    tail = t**(P+1)/math.factorial(P+1)/(1-t/(P+2))
    return ans, count, tail


def reduced(V, lam, tau, h, r, order):
    H = (V[:,r:]/lam[r:])@V[:,r:].T
    D = 1+1j*h*np.diag(H)
    R = np.exp(1j*np.pi/4)*np.sqrt(2*h)*(V[:,:r]/np.sqrt(lam[:r]))
    q,w = roots_hermite(order)
    w /= np.sqrt(np.pi)
    # iterate first axis, vectorize remaining r-1
    if r == 1:
        grid = np.empty((1,0)); ww = np.ones(1)
    else:
        grid = np.array(list(itertools.product(q,repeat=r-1)))
        ww = np.prod(np.array(list(itertools.product(w,repeat=r-1))),axis=1)
    ans = np.zeros(3,dtype=complex)
    for a,wa in zip(q,w):
        z = np.column_stack((np.full(len(grid),a),grid))@R.T
        frozen = np.exp(1j*tau*np.exp(-z*z/2).sum(axis=1))
        W = np.exp(-z*z/(2*D))/np.sqrt(D)
        F = np.exp(1j*tau*W.sum(axis=1))
        wp = -z/D*W
        pair = np.sum((wp@H)*wp,axis=1)
        corr = F*(1+1j*h*(1j*tau)**2*pair/2)
        ans += wa*np.array([ww@frozen,ww@F,ww@corr])
    return ans, h*h*tau*tau*np.exp(len(lam)*abs(tau))*np.sum(H*H)/4


def main():
    print(platform.python_version(),np.__version__,scipy.__version__)
    T,U,h=2.,.05,.1
    for N in (3,5,8):
        V,lam,G=covariance(N,T); tau=U*T/N/h
        start=time.perf_counter(); ref,terms,tail=reference(G,tau,h); secs=time.perf_counter()-start
        print(f'N={N} exact-series Z={ref:.13g} terms={terms} tail<={tail:.3g} seconds={secs:.4f}')
        for r in range(1,min(N-1,3)+1):
            a,_=reduced(V,lam,tau,h,r,24)
            start=time.perf_counter(); b,bound=reduced(V,lam,tau,h,r,40); elapsed=time.perf_counter()-start
            c,_=reduced(V,lam,tau,h,r,56)
            print(f' r={r} frozen,PA,corr errors={np.abs(c-ref)} bound={bound:.7g} delta24-40={max(abs(a-b)):.2g} delta40-56={max(abs(b-c)):.2g} nodes={40**r} seconds={elapsed:.4f}')
    # Actual large-N reduction only; no reference at N=160.
    V,lam,G=covariance(160,T); tau=U*T/160/h
    for r in (1,2,3):
        start=time.perf_counter(); b,bound=reduced(V,lam,tau,h,r,40); elapsed=time.perf_counter()-start
        c,_=reduced(V,lam,tau,h,r,56)
        print(f'N=160 r={r} Zcorr={c[2]:.13g} bound={bound:.7g} delta40-56={max(abs(b-c)):.2g} nodes={40**r} seconds={elapsed:.4f}')


if __name__=='__main__': main()
