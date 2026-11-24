import numpy as np
from scipy.special import gamma, roots_hermitenorm
from scipy.integrate import quad_vec, quad
from scipy.optimize import brentq
import itertools,sys

def run(N,h,kappa=0,ref_order=0):
    T=brentq(lambda T: N/2*(4*N/T*np.sin(np.pi/(2*N))**2-2*T/N)+kappa*h**(2/3),.1,3)
    d=T/N; tau=d/h; a=-kappa*h**(2/3)
    j=np.arange(1,N); k=np.arange(2,N)
    v=np.sin(np.pi*j/N); V=np.sqrt(2/N)*np.sin(np.pi*np.outer(j,k)/N)
    lam=4/d*np.sin(np.pi*k/(2*N))**2-2*d
    G=(V/lam)@V.T; C=h*G; diag=np.diag(C); sig=np.sqrt(diag)
    damp=lambda q: T/(8*h)*( (q*q/4)**3 if q*q<=2 else 3*q*q/16-1/4 )
    def w(r):
        t=r*r/2
        re=np.where(t<.01,-t**3/6+t**5/120-t**7/5040,np.sin(t)-t)
        return tau*(re-2j*np.sin(t/2)**2)
    def moment(k,p): return sig**k*(2**(k*p/2)*gamma((k*p+1)/2)/np.sqrt(np.pi))**(1/p)
    def norms(q,p):
        R=np.abs(q*v)
        A=tau*np.sum(R**5*moment(1,p)/8)
        B=tau*np.sum(3*R**2*moment(2,p)/4)
        D0=A+B; D=D0+tau*np.sum(R*moment(3,p)/2+moment(4,p)/8)
        T3=tau/6*np.sum((3*R+R**3)*moment(3,p)+(3+3*R**2)*moment(4,p)+3*R*moment(5,p)+moment(6,p))
        return D,D0,T3
    def bound(q):
        R=np.abs(q*v); D2,D02,T2=norms(q,2); D3,_,_=norms(q,3)
        U1=tau/24*np.sum((3+6*R**2+R**4)*moment(4,1)+(12*R+4*R**3)*moment(5,1)+(6+6*R**2)*moment(6,1)+4*R*moment(7,1)+moment(8,1))
        L=tau*np.sum(R**5*moment(1,1)/8+3*R**2*moment(2,1)/4+moment(4,1)/8)
        return np.array([L,U1+T2*(D2+D02)/2+D3**3/6])*np.exp(-a*q*q/(2*h)-damp(q))
    if ref_order:
        x,wt=roots_hermitenorm(ref_order); mesh=np.array(list(itertools.product(x,repeat=N-2)))
        weights=np.prod(np.array(list(itertools.product(wt/np.sqrt(2*np.pi),repeat=N-2))),axis=1)
        u=mesh@(V*np.sqrt(h/lam)).T
    def f(s):
        q=s*h**.25; r=q*v; t=r*r/2; E=np.exp(-1j*t)
        wp=tau*r*(E-1); wpp=tau*(E-1-1j*r*r*E); wppp=tau*E*(-3j*r-r**3)
        ell=wp.imag; al=wp.real; b=wpp/2; c=wppp/6
        mu=1j*(C@ell); fac=np.exp(-a*q*q/(2*h)+np.sum(w(r))-.5*ell@C@ell)
        EB=b@(mu**2+diag); EA=al@mu; cubic=c@(mu**3+3*mu*diag)
        delta2=(EA+EB)**2+al@C@al+4*(al@C)@(b*mu)+2*b@(C*C)@b+4*(b*mu)@C@(b*mu)
        prev=fac*(1+1j*b.imag@(mu**2+diag)); new=fac*(1+EA+EB+cubic+.5*delta2)
        out=[prev,new]
        if ref_order: out.append(np.exp(-a*q*q/(2*h))*(weights@np.exp(np.sum(w(r+u),axis=1))))
        return np.array(out)
    qmax=brentq(lambda q:damp(q)+a*q*q/(2*h)-65,.0001,10); smax=qmax/h**.25
    vals=2*quad_vec(f,0,smax,epsabs=1e-10,epsrel=1e-10,limit=10000)[0]
    bounds=2*quad_vec(lambda s:bound(s*h**.25),0,np.inf,epsabs=1e-10)[0]
    print(N,h,kappa,ref_order,'bounds',bounds,'values',vals,'errors',abs(vals[:2]-vals[2]) if ref_order else None,flush=True)
if __name__=='__main__':
    print(sys.version,np.__version__)
    for N in (3,4):
        for h in (1e-3,1e-4,1e-5): run(N,h,0,20)
    for h in (1e-4,1e-5,1e-6,1e-7): run(80,h)
    run(4,1e-4,.05,20); run(4,1e-4,.05,28)
