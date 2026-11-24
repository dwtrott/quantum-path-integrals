import numpy as np
from scipy.integrate import quad_vec, quad
from scipy.special import roots_hermitenorm
import sys, scipy
# a = first nonlinear Hessian coefficient in unnormalized sine coordinate.
def setup(N,h,kappa=0.,U=1.,m=1.,omega=1.):
    a=-kappa*h**(2/3); nu2=omega**2+U/m
    T=(-a+np.sqrt(a*a+4*m*m*nu2*N*N*np.sin(np.pi/(2*N))**2))/(m*nu2)
    d=T/N; j=np.arange(1,N); k=np.arange(2,N); v=np.sin(np.pi*j/N)
    E=np.sqrt(2/N)*np.sin(np.pi*j[:,None]*k[None,:]/N)
    lam=4*m/d*np.sin(k*np.pi/(2*N))**2-d*m*nu2
    assert np.min(lam)>0 and 4*m/d*np.sin(np.pi/(2*N))**2-d*m*omega**2>0
    L=E/np.sqrt(lam); Gdiag=np.sum(L*L,axis=1)
    return a,T,d,v,L,Gdiag

def run(N,h,kappa=0.,ref=False):
    U=1.; a,T,d,v,L,G=setup(N,h,kappa,U)
    tau=U*d/h; scale=h**.25
    def w(r):
        x=r*r/2
        return 1j*tau*np.where(np.abs(x)<.001,(-x*x/2+1j*x**3/6+x**4/24-1j*x**5/120-x**6/720),np.expm1(-1j*x)+1j*x).sum(axis=-1)
    def approx(q):
        r=q*v; ell=-tau*r*np.sin(r*r/2)
        return np.exp(-a*q*q/(2*h)+w(r)-h*np.sum((ell@L)**2)/2)
    def corrected(q):
        r=q*v; ell=-tau*r*np.sin(r*r/2)
        b=-tau*(np.sin(r*r/2)+r*r*np.cos(r*r/2))
        z=ell@L; gl=L@z
        return approx(q)*(1+.5j*np.sum(b*(h*G-h*h*gl*gl)))
    def frozen(q):return np.exp(-a*q*q/(2*h)+w(q*v))
    def damping(q):
        t=q*q/4
        return U*T/(8*h)*(t**3 if t<=.5 else .75*t-.25)
    C1=U*d*np.sum(v**5*np.sqrt(2*G/np.pi))/(8*np.sqrt(h)); C2=.75*U*d*np.sum(v*v*G); C0=.375*U*d*h*np.sum(G*G)
    def bfun(q):return np.exp(-a*q*q/(2*h)-damping(q))*(C1*q**5+C2*q*q+C0)
    # Bound integral evaluated approximately, not an interval-certified number.
    B=2*quad(lambda s:bfun(scale*s),0,np.inf,epsabs=1e-10,limit=300)[0]
    R=scale
    while damping(R)+a*R*R/(2*h)<42:R*=1.2
    J=2*quad_vec(lambda s:approx(scale*s),0,R/scale,epsabs=2e-9,epsrel=2e-9,limit=20000)[0]
    F=2*quad_vec(lambda s:frozen(scale*s),0,R/scale,epsabs=2e-9,epsrel=2e-9,limit=20000)[0]
    P=2*quad_vec(lambda s:corrected(scale*s),0,R/scale,epsabs=2e-9,epsrel=2e-9,limit=20000)[0]
    out=dict(N=N,h=h,kappa=kappa,T=T,a=a,J=J,F=F,P=P,B=B,cutoff=R,coeff=[C1,C2,C0])
    if ref:
        for order in (16,28):
            xx,ww=roots_hermitenorm(order); ww=ww/np.sqrt(2*np.pi); n=N-2
            mesh=np.meshgrid(*([xx]*n),indexing='ij'); weights=np.meshgrid(*([ww]*n),indexing='ij')
            us=np.stack([z.ravel() for z in mesh],axis=1)@L.T*np.sqrt(h)
            weights=np.prod(np.stack(weights,axis=-1),axis=-1).ravel()
            def fun(s):
                q=scale*s
                return np.exp(-a*q*q/(2*h))*np.sum(weights*np.exp(w(q*v+us)))
            I=2*quad_vec(fun,0,R/scale,epsabs=2e-9,epsrel=2e-9,limit=20000)[0]
            out[f'I{order}']=I
        out['error']=abs(J-I);out['frozen_error']=abs(F-I);out['corrected_error']=abs(P-I)
    print(out,flush=True)
if __name__=='__main__':
    print(sys.version, np.__version__, scipy.__version__,flush=True)
    for N in (3,4,80):
        for h in (1e-2,1e-3,1e-4,1e-5):run(N,h,0.,N<=4)
    for N in (3,80):run(N,1e-4,.05,N<=4)
