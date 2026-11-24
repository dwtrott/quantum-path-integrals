import numpy as np, scipy, sys
from scipy.special import roots_hermitenorm
from scipy.integrate import quad_vec
from scipy.optimize import brentq
from itertools import product

# Exact finite-dimensional model, m=omega=U=1. No certification.
def trial(N,h,k,order):
    j=np.arange(1,N); v=np.sin(np.pi*j/N)
    a=-k*np.sqrt(h)
    T=brentq(lambda T: N*N/T*(np.cos(2*np.arcsin(T*np.sqrt(2)/(2*N)))-np.cos(np.pi/N))-a, .8,3.0,xtol=1e-14)
    d=T/N; tau=d/h
    kk=np.arange(2,N); S=np.sqrt(2/N)*np.sin(np.pi*np.outer(j,kk)/N)
    lam=4/d*np.sin(np.pi*kk/(2*N))**2-2*d
    G=(S/lam)@S.T; C=h*G
    x,wt=roots_hermitenorm(order); wt/=np.sqrt(2*np.pi)
    idx=np.array(list(product(range(order),repeat=N-2)))
    us=x[idx]@((S*np.sqrt(h/lam)).T); ws=np.prod(wt[idx],axis=1)
    # Stable expm1 minus linear term via series at small |z|.
    def w(s):
        z=-.5j*s*s
        out=np.expm1(z)-z
        zz=z; p=z*z/2; ser=p.copy()
        for l in range(3,10): p=p*z/l; ser+=p
        return 1j*tau*np.where(abs(z)<.01,ser,out)
    scale=h**.25; rot=np.exp(-1j*np.pi/8)
    def integrand(s,angle,ref):
        q=scale*s*angle; r=q*v
        base=-a*q*q/(2*h)
        if ref: return angle*np.dot(ws,np.exp(base+np.sum(w(r+us),axis=1)))
        z=-.5j*r*r; wp=tau*r*np.expm1(z); wpp=tau*(np.expm1(z)-1j*r*r*np.exp(z))
        frozen=angle*np.exp(base+w(r).sum())
        return frozen*np.array([1,1+.5*(wpp@np.diag(C)+wp@C@wp)])
    # Real reference on large finite interval. Rotated approximation on Q=1/2.
    # Independent real reference, no cutoff in Gaussian expectation other than quadrature.
    R=2/scale
    ref,er=quad_vec(lambda s: integrand(s,1,True),0,R,epsabs=2e-11,epsrel=2e-11,limit=2000)
    Q=.5/scale
    pred,ep=quad_vec(lambda s: integrand(s,rot,False),0,Q,epsabs=2e-12,epsrel=2e-12,limit=1000)
    return 2*ref,2*pred,er,ep
if __name__=='__main__':
    print(sys.version, np.__version__,scipy.__version__)
    for N in [3,4]:
        for h in [1e-4,1e-5,1e-6]:
            for k in [-.5,.5]:
                r,p,er,ep=trial(N,h,k,18); r2,_,_,_=trial(N,h,k,28)
                print(N,h,k,'ref',r,'refdiff',abs(r-r2),'errors',abs(p-r2),'quad',er,ep,flush=True)
