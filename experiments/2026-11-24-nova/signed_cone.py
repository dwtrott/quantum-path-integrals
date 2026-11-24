"""Signed directional curvature with a mean cone; run OPENBLAS_NUM_THREADS=1 python signed_cone.py."""
import sys,time
import numpy as np
import scipy,mpmath
from mpmath import iv
from scipy.linalg import solve
from structured_audit import setup
from structured_middle_obstruction import root_ball,I
iv.dps=45

def certify(N,L=64):
    start=time.perf_counter();_,_,_,K,D,_,z,_,_=setup(N);n=N+1;d=I(2)/N;rho=I('1e-8')
    m=[d/2 if j in (0,N) else d for j in range(n)]
    zz=[v+iv.mpc(iv.mpf([-1,1])*rho,iv.mpf([-1,1])*rho) for v in root_ball(N,z,K,D,m,d,rho)]
    s=I('.8291');a=I('-1.6247');r=I('.49');u=1+1j*s;b=1j*a;gamma=I(2);cap=I('.05');delta=I('.04')
    assert (s-I('.02'))*r<I('.4') and iv.sqrt(2)*cap<I('.4')
    # Rectangles for the radial coordinate x and scaled mean y. At radial panel k,
    # x in [-kr/L,kr/L], y in [-k cap/L,k cap/L]. Bound by a cumulative x sweep
    # with the full mean cap (looser, but O(NL)); split mean interval into 8 parts.
    lower=[];cross=[];mean=[]
    for v in zz:
        mn=None;mx=None;cn=None;acc=I(0);acb=I(0);accm=I(0)
        for k in range(1,L+1):
            l=r*(k-1)/L;h=r*k/L
            for xx in [iv.mpf([l.a,h.b]),iv.mpf([-h.b,-l.a])]:
                for j in range(8):
                    yl=cap*(j-4)/4;yh=cap*(j-3)/4
                    w=v+u*xx+b*iv.mpf([yl.a,yh.b]);g=3*(1-w*w)*iv.exp(-w*w/2)
                    A=(u*u*g).imag;B=(u*b*g).imag;C=(b*b*g).imag
                    mn=A.a if mn is None else min(mn,A.a)
                    mx=abs(B).b if mx is None else max(mx,abs(B).b)
                    cn=C.a if cn is None else min(cn,C.a)
            wt=I(2*(L-k)+1)/(2*L*L)
            acc+=wt*mn;acb+=wt*mx;accm+=wt*cn
        lower.append(acc.a);cross.append(acb.b);mean.append(accm.a)
    off=-(s-delta)/d
    def test(lam):
        diag=[(s-delta)*(1 if j in (0,N) else 2)/d+((1-s*s)/2 if j in (0,N) else 0)+m[j]*(lower[j]-lam) for j in range(n)]
        piv=[diag[0]];f=[m[0]]
        for j in range(1,n):
            if not (piv[-1]>0 or piv[-1]<0):return False,None
            f.append(m[j]-off/piv[-1]*f[-1]);piv.append(diag[j]-off**2/piv[-1])
        if not all(p>0 or p<0 for p in piv):return False,None
        schur=-1/gamma-sum(x*x/p for x,p in zip(f,piv))
        return sum(p<0 for p in piv)==1 and schur>0,schur
    # Exact discrete trace: |v_0+v_N-2*mu| <= sqrt((2/3)*(1-1/N**2)*E).
    G=iv.sqrt(sum(m[j]*cross[j]**2 for j in range(n)))
    C=4*gamma+a*a+2*s*a+(s*a)**2*(I(2)/3)*(1-I(1)/N**2)/(4*delta)-sum(m[j]*mean[j] for j in range(n))
    for lam0 in ['.4','.5','.6','.7','.8']:
        lam=I(lam0);ok,schur=test(lam)
        for eta0 in ['.04','.05','.06','.07']:
            eta=I(eta0);assert eta*iv.sqrt(2)*r<cap
            margin=lam-2*G*eta-max(C.b,I(0))*eta**2
            print('N',N,'L',L,'lambda',lam0,'verified',ok,'eta',eta0,'margin lower',float(margin.a),'G upper',float(G.b),flush=True)
            if lam0=='.5' and eta0=='.07':assert ok and margin>I('.15')
    print('seconds',time.perf_counter()-start,flush=True)
if __name__=='__main__':
    print('versions',np.__version__,scipy.__version__,mpmath.__version__,flush=True)
    for N in map(int,sys.argv[1:] or [16,32,64,128]):certify(N)
