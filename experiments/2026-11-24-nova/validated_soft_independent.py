import numpy as np, scipy, sys
from scipy.integrate import quad_vec
for N,eta in [(4,0),(80,.001),(80,-.001)]:
    h=1e-5; th=(np.pi+eta)/N; d=np.sqrt(2)*np.sin(th/2); a=-2*N/d*np.sin((np.pi+eta/2)/N)*np.sin(eta/(2*N))
    j=np.arange(1,N); k=np.arange(2,N); v=np.sin(j*np.pi/N)
    S=np.sqrt(2/N)*np.sin(np.pi*np.outer(j,k)/N); C=h*(S/((2/d)*(np.cos(th)-np.cos(k*np.pi/N))))@S.T
    rot=np.exp(-1j*np.pi/8); scale=h**.25
    def f(x):
        q=scale*rot*x; r=q*v; z=-.5j*r*r; E=np.exp(z); tau=d/h
        w=tau*1j*np.sum(np.expm1(z)-z); wp=tau*r*np.expm1(z); wpp=tau*(np.expm1(z)-1j*r*r*E)
        return rot*np.exp(-a*q*q/(2*h)+w)*(1+(wpp@np.diag(C)+wp@C@wp)/2)
    value,err=quad_vec(f,0,.5/scale,epsabs=1e-12,epsrel=1e-12)
    print(N,eta,2*value,2*err)
print(sys.version,np.__version__,scipy.__version__)
