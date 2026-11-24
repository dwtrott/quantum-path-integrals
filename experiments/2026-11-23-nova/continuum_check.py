"""Conventional continuum certificate, numerical checks of p=1,2 only.
Run python continuum_check.py; scipy quadrature and ordinary float, not intervals.
"""
import sys, math, numpy as np, scipy
from scipy.integrate import quad
from math import exp, sin, sqrt, asin, pi
h=.1; U=.05; m=omega=1.
L=[1.,sqrt(2/pi),1.,2*sqrt(2/pi),3.]

def effective(T,N):
    d=T/N; theta=2*asin(omega*d/2)
    return theta/d, m*sin(theta)/theta

def bridge(s,t,T,nu,M):
    return np.sin(nu*np.minimum(s,t))*np.sin(nu*(T-np.maximum(s,t)))/(M*nu*sin(nu*T))

def path(s,T,nu,a,b):
    return (a*np.sin(nu*(T-s))+b*np.sin(nu*s))/sin(nu*T)

def moments(t):
    a,b,c,d=[t*x for x in L[1:]]
    return exp(t)*np.array([1,a,b+a*a,c+3*a*b+a**3,d+4*a*c+3*b*b+6*a*a*b+a**4])

def fixed_p_moments(p):
    # Moments of a sum of p independent |standard normal| variables.
    v=np.array([1.,0,0,0,0])
    for _ in range(p):
        v=np.array([sum(math.comb(j,k)*v[k]*L[j-k] for k in range(j+1)) for j in range(5)])
    return v

def budgets(T,N,a,b,B,tuned=False):
    d=T/N; nu,M=(omega,m) if tuned else effective(T,N); E=abs(a)+abs(b)
    eps=nu-omega; s0=abs(sin(omega*T)); s=s0-T*eps
    if s<=0: return float('inf'),float('inf')
    se=abs(sin(nu*T)); X1=E*nu/se; X2=E*nu**2/se
    G1=1/(M*se); G2=nu/(M*se)
    _,b1,b2,b3,b4=B
    quaderr=d*d/8*(X2*b1+(X1*X1+2*h*G2+h/(M*T))*b2+2*X1*h*G1*b3+h*h*G1*G1*b4)
    dx=E*T*eps*(1/s+1/s**2)
    dg=abs(1/(M*nu)-1/(m*omega))/s+T*eps*(1/s+1/s**2)/(m*omega)
    paramerr=dx*b1+h*dg*b2/2
    return quaderr,paramerr

def f1(s,T,nu,M,a,b):
    x=path(s,T,nu,a,b); A=1+1j*h*bridge(s,s,T,nu,M)
    return np.exp(-x*x/(2*A))/np.sqrt(A)

def f2(s,t,T,nu,M,a,b):
    x=path(s,T,nu,a,b); y=path(t,T,nu,a,b)
    A=1+1j*h*bridge(s,s,T,nu,M); D=1+1j*h*bridge(t,t,T,nu,M); C=1j*h*bridge(s,t,T,nu,M)
    det=A*D-C*C
    return np.exp(-(D*x*x-2*C*x*y+A*y*y)/(2*det))/np.sqrt(det)

def integrate(f,a,b):
    return complex(quad(lambda x:f(x).real,a,b,epsabs=1e-12,epsrel=1e-12)[0],quad(lambda x:f(x).imag,a,b,epsabs=1e-12,epsrel=1e-12)[0])

def ref(T,a,b):
    return [integrate(lambda s:f1(s,T,omega,m,a,b),0,T)/T,2*integrate(lambda t:integrate(lambda s:f2(s,t,T,omega,m,a,b),0,t),0,T)/T**2]

def grid(T,N,a,b,tuned=False):
    nu,M=(omega,m) if tuned else effective(T,N); x=np.linspace(0,T,N+1); w=np.ones(N+1)/N; w[[0,-1]]*=.5
    return [w@f1(x,T,nu,M,a,b),w@f2(x[:,None],x[None,:],T,nu,M,a,b)@w]

if __name__=='__main__':
    print(sys.version.split()[0],np.__version__,scipy.__version__)
    for T,a,b in [(2.,0.,0.),(3.,.2,-.1),(3.3,.2,-.1)]:
        z=ref(T,a,b); print('T, a, b, coefficient references',T,a,b,z)
        for tuned in (False,True):
            print('tuned',tuned,'N, |delta c1|, bound1, |delta c2|, bound2, full normalized bound')
            for N in (16,32,64,128,256,640,2560):
                zz=grid(T,N,a,b,tuned) if N<=256 else None
                b1=sum(budgets(T,N,a,b,fixed_p_moments(1),tuned)); b2=sum(budgets(T,N,a,b,fixed_p_moments(2),tuned))
                full=sum(budgets(T,N,a,b,moments(abs(U)*T/h),tuned))
                if zz is not None:
                    e1,e2=abs(zz[0]-z[0]),abs(zz[1]-z[1]); assert e1<=b1 and e2<=b2
                    print(N,*(f'{v:.8g}' for v in [e1,b1,e2,b2,full]))
                else: print(N,'-','-','-','-',f'{full:.8g}')
