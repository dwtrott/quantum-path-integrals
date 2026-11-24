import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import scipy
from scipy.integrate import quad_vec
import time, sys

# eta=N*(theta-pi/N), m=omega=U=1. d=2*sin(theta/2)/sqrt(2).
def generators(N,eta):
    alpha=np.pi/N; delta=eta/N; theta=alpha+delta
    d=np.sqrt(2)*np.sin(theta/2); j=np.arange(1,N,dtype=float); k=N-j
    v=np.sin(j*alpha); vk=np.sin(k*alpha)
    sinc=lambda t:np.sinc(t/np.pi)
    S=sinc(eta); s=sinc(delta/2); mid=alpha+delta/2
    K=1/(N*np.sin(theta)*S)
    # (sinc(eta)-sinc(delta/2))/delta; small eta series
    if abs(eta)<.5:
        ds=0.; fact=1.
        for r in range(1,10):
            fact*=(2*r)*(2*r+1)
            ds+=(-1)**r*N*eta**(2*r-1)*(1-(.5/N)**(2*r))/fact
    else:
        ds=(S-s)/delta
    D=(.5*np.cos(alpha+3*delta/4)*sinc(delta/4)*S+np.sin(mid)*ds)/(N*np.sin(mid)*s*np.sin(theta)*S)
    J=j*np.cos(j*(alpha+delta/2))*sinc(j*delta/2)
    L=k*np.cos(k*(alpha+delta/2))*sinc(k*delta/2)
    P=np.stack([d*v,-d*K*J],axis=1)
    Q=np.stack([D*vk-K*L,np.sin(k*theta)],axis=1)
    diag=np.sum(P*Q,axis=1)
    a=-2*N/d*np.sin(mid)*np.sin(delta/2)
    return d,a,v,P,Q,diag

def mv(P,Q,x):
    # symmetric, diagonal counted exactly once
    A=P*x[:,None]; B=Q*x[:,None]
    left=np.cumsum(A,axis=0); right=np.cumsum(B[::-1],axis=0)[::-1]
    return np.sum(Q*left+P*(right-B),axis=1)

def square_gen(P,Q):
    return (np.stack([P[:,0]**2,2*P[:,0]*P[:,1],P[:,1]**2],axis=1),
            np.stack([Q[:,0]**2,Q[:,0]*Q[:,1],Q[:,1]**2],axis=1))

def dense(N,d,eta):
    j=np.arange(1,N); k=np.arange(2,N)
    V=np.sqrt(2/N)*np.sin(np.pi*np.outer(j,k)/N)
    theta=(np.pi+eta)/N
    lam=4/d*np.sin((np.pi*k/N+theta)/2)*np.sin((np.pi*k/N-theta)/2)
    return (V/lam)@V.T

def evaluator(N,h,eta,fast):
    d,a,v,P,Q,diag=generators(N,eta); tau=d/h
    if fast:
        P=P*h; diag=diag*h; P2,Q2=square_gen(P,Q)
        M=lambda x:mv(P,Q,x); M2=lambda x:mv(P2,Q2,x)
    else:
        C=h*dense(N,d,eta); diag=np.diag(C); CC=C*C
        M=lambda x:C@x; M2=lambda x:CC@x
    def f(s):
        q=s*h**.25; r=q*v; t=r*r/2; E=np.exp(-1j*t)
        w=tau*(np.where(t<.01,-t**3/6+t**5/120-t**7/5040,np.sin(t)-t)-2j*np.sin(t/2)**2)
        wp=tau*r*(E-1); b=tau*(E-1-1j*r*r*E)/2; c=tau*E*(-3j*r-r**3)/6
        ell=wp.imag; al=wp.real; Cell=M(ell); mu=1j*Cell
        fac=np.exp(-a*q*q/(2*h)+np.sum(w)-.5*ell@Cell)
        bm=b*mu; EA=al@mu; EB=b@(mu**2+diag); cubic=c@(mu**3+3*mu*diag)
        delta2=(EA+EB)**2+al@M(al)+4*al@M(bm)+2*b@M2(b)+4*bm@M(bm)
        return fac*(1+EA+EB+cubic+.5*delta2)
    return f

def run():
    print(sys.version, np.__version__, scipy.__version__)
    rng=np.random.default_rng(37)
    for N in (3,4,80,640):
        for eta in (-.2,-1e-12,0.,1e-12,.2):
            d,a,v,P,Q,di=generators(N,eta); G=dense(N,d,eta); x=rng.normal(size=N-1)+1j*rng.normal(size=N-1)
            P2,Q2=square_gen(P,Q)
            err=max(np.max(abs(di-np.diag(G))),np.max(abs(mv(P,Q,x)-G@x))/max(1,np.max(abs(G@x))),np.max(abs(mv(P2,Q2,x)-(G*G)@x))/max(1,np.max(abs((G*G)@x))))
            assert err<1e-10,(N,eta,err)
            print('covariance',N,eta,err)
    # compare actual integrand, same quadrature on fixed interval, no full kernel reference
    for N in (4,80,640):
        h=1e-5; eta=0.; out=[]
        for fast in (False,True):
            t=time.perf_counter(); f=evaluator(N,h,eta,fast); setup=time.perf_counter()-t
            t=time.perf_counter(); val=2*quad_vec(f,0,20,epsabs=1e-8,epsrel=1e-8,limit=10000)[0]; elapsed=time.perf_counter()-t
            out.append(val); print('integral',N,fast,setup,elapsed,val,flush=True)
        if len(out)==2: assert abs(out[0]-out[1])<1e-7
    for N in (640,2560,10240,40960):
        for fast in ([False,True] if N<=2560 else [True]):
            t=time.perf_counter(); f=evaluator(N,1e-5,0.,fast); setup=time.perf_counter()-t
            t=time.perf_counter(); vals=[f(1.2) for _ in range(10)]; elapsed=(time.perf_counter()-t)/10
            print('sample',N,fast,setup,elapsed,vals[-1],flush=True)
    # negative result: subtracting singular inverse loses digits even for tiny N.
    N=80; eta=1e-12; d,a,v,P,Q,di=generators(N,eta); j=np.arange(1,N); th=(np.pi+eta)/N
    raw=(d*np.sin(np.minimum.outer(j,j)*th)*np.sin((N-np.maximum.outer(j,j))*th)/(np.sin(th)*np.sin(N*th)))-np.outer(v,v)/a
    print('naive subtraction max error',np.max(abs(raw-dense(N,d,eta))))
if __name__=='__main__': run()
