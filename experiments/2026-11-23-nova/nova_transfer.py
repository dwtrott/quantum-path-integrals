"""Nova: matched finite-N benchmark. Floating-point only; analytic errors exclude roundoff.
Run: OPENBLAS_NUM_THREADS=1 python nova_transfer.py
"""
import math, sys, time
import numpy as np
import scipy
from scipy.signal import fftconvolve
from scipy.fft import dst, idst

T,U,h,m,omega=2.,.05,.1,1.,1.

def params(N):
    d=T/N; n=N-1; k=np.arange(1,N)
    lam=m*d*(4/d**2*np.sin(np.pi*k/(2*N))**2-omega**2)
    # Conservative upper bound on largest Gaussian marginal variance.
    vmax=2*h/N*np.sum(1/lam)
    return d,lam,vmax

def transfer(N, eps=1e-8, L=18, use_fft=True):
    d,lam,vmax=params(N); n=N-1; tau=U*d/h; t=n*abs(tau)
    Lambda=4*m/(d*h)+L**2*vmax
    # Give alias and box each <= eps/4, leaving room for Taylor tail.
    a=math.log(32*n*math.exp(t)/eps)
    s=math.sqrt(2*math.pi**2/(Lambda*a))
    M=math.ceil(math.sqrt(2*vmax*math.log(16*n*math.exp(t)/eps))/s)
    v=math.exp(-a); logtheta=math.log1p(2*v/(1-v))
    alias=math.exp(t)*math.expm1(n*logtheta)
    R=t**(L+1)/math.factorial(L+1)/(1-t/(L+2))
    tay=R*(1+math.exp(n*logtheta))
    box=math.exp(t+n*logtheta)*2*n*math.exp(-(s*M)**2/(2*vmax))
    grid=s*np.arange(-M,M+1)
    c=m/(2*h*d); pref=math.sqrt(c/math.pi)
    weight=np.exp(d*m*omega**2*grid**2/(2*h)+1j*tau*np.exp(-1j*grid**2/2))
    f=pref*np.exp(-c*grid**2)*weight
    kernel=s*pref*np.exp(-c*(s*np.arange(-2*M,2*M+1))**2)
    for _ in range(1,n):
        f=weight*(fftconvolve(f,kernel,mode='same') if use_fft else np.convolve(f,kernel,mode='full')[2*M:4*M+1])
    num=s*np.dot(pref*np.exp(-c*grid**2),f)
    den=pref*np.exp(-.5*np.sum(np.log(lam*d/m)))
    return num/den,alias+tay+box,2*M+1,s, (alias,tay,box)

def pa(N,r=2,order=40):
    d,lam,_=params(N); n=N-1; tau=U*d/h; t=n*abs(tau)
    j=np.arange(1,N); v=np.sqrt(2/N)*np.sin(np.pi*j[:,None]*np.arange(1,r+1)[None,:]/N)
    a=np.zeros(N); a[r+1:]=1/lam[r:]
    diag=(np.sum(a)-np.fft.fft(a).real[1:N])/N
    D=1+1j*h*diag; zfactor=np.exp(1j*np.pi/4)*np.sqrt(2*h)*v/np.sqrt(lam[:r])
    xs,ws=np.polynomial.hermite.hermgauss(order)
    inv=np.zeros(n);inv[r:]=1/lam[r:]
    result=0j; prim=0j
    for x,wx in zip(xs,ws):
        q=np.column_stack((np.full(order,x),xs)); z=q@zfactor.T
        W=np.exp(-z*z/(2*D))/np.sqrt(D); wp=-z*W/D
        F=np.exp(1j*tau*W.sum(axis=1))
        wdst=dst(wp,type=1,norm='ortho',axis=1)
        P=1j*h*(1j*tau)**2*np.sum(wdst*wdst*inv,axis=1)/2
        result+=wx*np.dot(ws,F*(1+P)); prim+=wx*np.dot(ws,F)
    bound=h*h*tau*tau*np.exp(t)*np.sum(inv*inv)/4
    return result/np.pi,prim/np.pi,bound

if __name__=='__main__':
    print(sys.version.split()[0],np.__version__,scipy.__version__)
    for N in (8,40,160,640,2560):
        start=time.perf_counter(); z,b,K,s,parts=transfer(N); t1=time.perf_counter()-start
        start=time.perf_counter(); p,prim,bp=pa(N); t2=time.perf_counter()-start
        print(N, 'transfer',z,'bound',b,'grid',K,'seconds',t1, 'PA',p,'PA bound',bp,'actual diff',abs(z-p),'PA seconds',t2,flush=True)
        print('parts',parts,'primitive diff',abs(z-prim),flush=True)
        if N<=160:
            zz,bb,_,_,_=transfer(N,eps=1e-10,use_fft=False)
            print('direct convolution refined diff',abs(z-zz), 'refined bound',bb,flush=True)
