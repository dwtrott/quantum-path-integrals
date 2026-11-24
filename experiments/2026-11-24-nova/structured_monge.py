"""Exact integer resource DP on a bounded uniform lattice, not a continuum certificate.
Run OPENBLAS_NUM_THREADS=1 python structured_monge.py
"""
import time
from itertools import product
import numpy as np
import numba
from numba import njit
from mpmath import iv
import mpmath, scipy
from fractions import Fraction as F
from structured_audit import setup
from structured_middle_obstruction import root_ball,I
iv.dps=60
INF=np.int64(1<<60)
SCALE=10**8

@njit
def rowmin(a,e,finite,lo,hi,left,right,out):
    if left>right:return
    mid=(left+right)//2
    best=INF; besti=lo
    for i in range(lo,hi+1):
        j=finite[i]; val=a[j]+e[j,mid]
        if val<best:best=val;besti=i
    out[mid]=best
    rowmin(a,e,finite,lo,besti,left,mid-1,out)
    rowmin(a,e,finite,besti,hi,mid+1,right,out)

@njit
def solve_dp(site,edges,threshold,fast):
    n,J=site.shape; K=J//2; S=2*(n-1)*K; W=2*S+1
    old=np.full((2,W,J),INF,dtype=np.int64)
    for k in range(J):old[int(abs(k-K)>=threshold),S+k-K,k]=site[0,k]
    reach=K
    for j in range(1,n):
        wt=1 if j==n-1 else 2
        new=np.full_like(old,INF)
        for flag in range(2):
            for s in range(S-reach,S+reach+1):
                a=old[flag,s]; finite=np.where(a<INF)[0]
                if len(finite)==0:continue
                out=np.full(J,INF,dtype=np.int64)
                if fast:rowmin(a,edges[j-1],finite,0,len(finite)-1,0,J-1,out)
                else:
                    for k in range(J):
                        for i in finite:out[k]=min(out[k],a[i]+edges[j-1,i,k])
                for k in range(J):
                    f=flag or abs(k-K)>=threshold
                    new[f,s+wt*(k-K),k]=min(new[f,s+wt*(k-K),k],out[k]+site[j,k])
        old=new;reach+=wt*K
    lo=np.empty(W,np.int64);loc=np.empty(W,np.int64)
    for s in range(W):lo[s]=np.min(old[1,s]);loc[s]=np.min(old[0,s])
    return lo,loc

def costs(N,K,step):
    d,m,E,A,D,t,z,phi,grad=setup(N); n=N+1; di=I(2)/N; rho=I('1e-8'); mi=[di/2 if j in (0,N) else di for j in range(n)]
    zc=root_ball(N,z,A,D,mi,di,rho)
    zz=[w+iv.mpc(iv.mpf([-1,1])*rho,iv.mpf([-1,1])*rho) for w in zc]
    X=[k*step for k in range(-K,K+1)]; J=len(X); sig=F('0.8291'); B=F('0.4')
    t=[F('.02')+F('1.4')*F(j,N)*(1-F(j,N)) for j in range(n)]
    Y=[[tj*x+max(-B,min(B,(sig-tj)*x)) for x in X] for tj in t]
    edges=np.empty((N,J,J),np.int64)
    for j in range(N):
        for i,x in enumerate(X):
            for k,x1 in enumerate(X):
                v=(x1-x)*(Y[j+1][k]-Y[j][i])*F(N,2)*SCALE
                edges[j,i,k]=v.numerator//v.denominator
        # Exact integer Monge check for all adjacent quadrangles; implies all quadrangles.
        e=edges[j];assert np.all(e[:-1,:-1]+e[1:,1:]<=e[:-1,1:]+e[1:,:-1])
    site=np.empty((n,J),np.int64); upper=np.empty_like(site)
    W=lambda z:iv.exp(-z*z/2)
    rat=lambda x:I(x.numerator)/I(x.denominator)
    for j in range(n):
        for k,x in enumerate(X):
            xx=rat(x);yy=rat(Y[j][k]);w=iv.mpc(xx,yy)
            g=3*mi[j]*(W(zz[j])-W(zz[j]+w)-zz[j]*W(zz[j])*w).imag
            if j in (0,N):g+=(xx**2-yy**2)/2
            # Convert an approximate endpoint, then validate direction using interval comparisons.
            l=int(np.floor(float((g*SCALE).a)))-1;u=int(np.ceil(float((g*SCALE).b)))+1
            assert I(l)<=g*SCALE and I(u)>=g*SCALE
            site[j,k]=l;upper[j,k]=u
    # Ensure all possible finite sums far below the sentinel and int64 overflow.
    assert (n+N)*max(np.max(abs(site)),np.max(abs(upper)),np.max(abs(edges)))<INF//4
    return site,edges,upper

def run(N,K,step,threshold,compare):
    start=time.perf_counter();site,edges,upper=costs(N,K,step);prep=time.perf_counter()-start
    start=time.perf_counter();lo,local=solve_dp(site,edges,threshold,True);fast=time.perf_counter()-start
    naive=None
    if compare:
        start=time.perf_counter();r,rl=solve_dp(site,edges,threshold,False);naive=time.perf_counter()-start
        assert np.array_equal(r,lo) and np.array_equal(rl,local)
    S=(len(lo)-1)//2
    # Upper bound from the same minimization with enclosing costs. Edge floor +1 is an upper bound.
    hi,_=solve_dp(upper,edges+1,threshold,True)
    assert np.all(lo<=hi)
    print('N',N,'K',K,'step',str(step),'outside grid index',threshold,'setup_s',prep,'fast_s',fast,'naive_s',naive,'exact-zero-mean bounds',lo[S]/SCALE,hi[S]/SCALE,'mean tolerance 4 lower',min(lo[S-4:S+5])/SCALE,'unconstrained lower',min(lo)/SCALE,flush=True)

if __name__=='__main__':
    print('versions',np.__version__,scipy.__version__,mpmath.__version__,numba.__version__,flush=True)
    # Exhaustive control for indexing, flag, sentinel and negative site costs.
    rng=np.random.default_rng(41);a=rng.integers(-50,50,(4,5),dtype=np.int64)
    x=np.arange(5,dtype=np.int64)-2;e=np.array([((x[:,None]-x[None,:])**2)*(j+1) for j in range(3)])
    brute=np.full((2,25),INF,np.int64)
    for v in product(range(5),repeat=4):
        s=v[0]+2*v[1]+2*v[2]+v[3];f=int(any(abs(k-2)>=2 for k in v))
        cost=sum(a[j,k] for j,k in enumerate(v))+sum(e[j,v[j],v[j+1]] for j in range(3))
        brute[f,s]=min(brute[f,s],cost)
    for fast in (True,False):
        lo,loc=solve_dp(a,e,2,fast);assert np.array_equal(lo,brute[1]) and np.array_equal(loc,brute[0])
    print('exhaustive control: 625 paths passed',flush=True)
    for args in [(16,10,F('0.1'),3,True),(16,20,F('0.1'),3,True),(16,40,F('0.05'),6,True),(32,40,F('0.05'),6,True),(32,80,F('0.025'),12,True)]:run(*args)
