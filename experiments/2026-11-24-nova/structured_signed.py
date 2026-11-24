"""Signed integrated Hessian and mean-cone certificate. OPENBLAS_NUM_THREADS=1."""
import sys,time
import numpy as np
import scipy,mpmath
from mpmath import iv
from structured_audit import setup
from structured_middle_obstruction import root_ball,I
iv.dps=45

def certify(N,L=128):
    start=time.perf_counter();_,_,_,K,D,_,z,_,_=setup(N);n=N+1;d=I(2)/N;rho=I('1e-8')
    m=[d/2 if j in (0,N) else d for j in range(n)]
    zz=[v+iv.mpc(iv.mpf([-1,1])*rho,iv.mpf([-1,1])*rho) for v in root_ball(N,z,K,D,m,d,rho)]
    sigma=I('.8291');a=I('-1.6247');r=I('.49');u=1+1j*sigma;delta=I('.08');gamma=I(2)
    assert (sigma-I('.02'))*r<I('.4') and max(abs(v.imag).b for v in zz)<I('.5')
    # For t in [(k-1)/L,k/L], tx belongs to [-kr/L,kr/L].
    # Running minima of signed Hessian along this expanding interval.
    lower=[]
    for v in zz:
        mn=None;acc=I(0)
        for k in range(1,L+1):
            l=r*(k-1)/L;h=r*k/L
            for xx in [iv.mpf([l.a,h.b]),iv.mpf([-h.b,-l.a])]:
                w=v+u*xx;f=(u*u*3*(1-w*w)*iv.exp(-w*w/2)).imag
                mn=f.a if mn is None else min(mn,f.a)
            wt=I(2*(L-k)+1)/(2*L*L) # integral of (1-t) over this t panel
            acc+=wt*mn
        lower.append(acc.a)
    # Q = (sigma-delta)K + endpoint quadratic + diag(m*lower).
    # Certify Q-lambda M+gamma*m*m^T via bordered tridiagonal inertia.
    off=-(sigma-delta)/d
    def test(lam):
        diag=[(sigma-delta)*(1 if j in (0,N) else 2)/d+((1-sigma*sigma)/2 if j in (0,N) else 0)+m[j]*(lower[j]-lam) for j in range(n)]
        piv=[diag[0]];f=[m[0]]
        for j in range(1,n):
            if not (piv[-1]>0 or piv[-1]<0):return False,None
            f.append(m[j]-off/piv[-1]*f[-1]);piv.append(diag[j]-off**2/piv[-1])
        if not all(p>0 or p<0 for p in piv):return False,None
        schur=-1/gamma-sum(x*x/p for x,p in zip(f,piv))
        return sum(p<0 for p in piv)==1 and schur>0,schur
    for lam in ['.7','.75']:
        ok,schur=test(I(lam));print('N',N,'L',L,'lambda',lam,'verified',ok,'Schur lower',None if schur is None else float(schur.a),flush=True)
        if lam=='.7':assert ok
    # Cone: |mu| <= eta*||v||_M, r=.49 implies |mu|<.012.
    eta=I('.016');lam=I('.7');G=I(11)
    assert eta*iv.sqrt(2)*r<I('.012') and iv.sqrt(2)*eta*iv.sqrt(2)*r<I('.4')
    assert I('.5')+sigma*r+abs(a)*I('.012')<1
    assert 3*(2*iv.exp(-I('.5'))+1)*iv.exp(I('.5'))<G
    linear=G*abs(a)*iv.sqrt(2)*abs(u)
    quadratic=4*gamma+G*a*a+2*sigma*abs(a)+a*a+2*sigma*sigma*a*a/delta
    cone=lam-linear*eta-quadratic*eta*eta
    assert cone>I('.1')
    print('N',N,'cone lower',float(cone.a),'seconds',time.perf_counter()-start,flush=True)

if __name__=='__main__':
    print('versions',np.__version__,scipy.__version__,mpmath.__version__,flush=True)
    for N in map(int,sys.argv[1:] or [16,32,64,128]):certify(N)
