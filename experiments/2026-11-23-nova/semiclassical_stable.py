# Nova, 2026-11-23. Run: OPENBLAS_NUM_THREADS=1 python semiclassical_stable.py
import sys,time
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import roots_hermitenorm
import scipy

def green_diag(N,T,U,m=1.,omega=1.):
 d=T/N;nu=np.sqrt(omega**2+U/m);theta=2*np.arcsin(d*nu/2)
 assert N*theta<np.pi
 j=np.arange(1,N);return d*np.sin(j*theta)*np.sin((N-j)*theta)/(m*np.sin(theta)*np.sin(N*theta))

def approximations(N,T,U,h):
 d=T/N;G=green_diag(N,T,U);s2=np.sum(G**2)*d;s3=np.sum(G**3)*d
 c=3*U*h*s2/8
 b2=h*h*(15*U*s3/48+105*U**2*s2*s2/128)
 return 1-1j*c,c,b2,np.max(G)

def second_correction(N,T,U,h):
 d=T/N;theta=2*np.arcsin(d*np.sqrt(1+U)/2);j=np.arange(1,N)
 left=np.sin(j*theta);right=np.sin((N-j)*theta);c0=d/(np.sin(theta)*np.sin(N*theta));G=c0*left*right
 S2=d*np.sum(G**2);S3=d*np.sum(G**3);S4=d*np.sum(G**4)
 def off(a,b):return 2*np.dot(np.cumsum(a)[:-1],b[1:])
 Q2=d*d*(np.sum(G**4)+c0*c0*off(G*left**2,G*right**2))
 Q4=d*d*(np.sum(G**4)+c0**4*off(left**4,right**4))
 EA2=U*U*h*h/64*(9*S2*S2+72*Q2+24*Q4)
 p=1-1j*3*U*h*S2/8-15*U*h*h*S3/48-EA2/2
 b3=h**3*(105*U*S4/384+945*U*U*S2*S3/384+10395*U**3*S2**3/3072)
 return p,b3

def transfer(N,T,U,h,grid=301,order=32,width=10):
 d=T/N;theta=2*np.arcsin(d*np.sqrt(1+U)/2);G=green_diag(N,T,U)
 x=np.linspace(-width,width,grid);nodes,weights=roots_hermitenorm(order);weights/=np.sqrt(2*np.pi)
 def weight(j,y):
  v=h*G[j]*y*y/2
  return np.exp(1j*U*d/h*(np.expm1(-1j*v)+1j*v))
 v=weight(N-2,x)
 for j in range(N-3,-1,-1):
  k=j+1
  corr=np.sqrt(np.sin(k*theta)*np.sin((N-k-1)*theta)/(np.sin((k+1)*theta)*np.sin((N-k)*theta)))
  y=corr*x[:,None]+np.sqrt(1-corr*corr)*nodes[None,:]
  vals=CubicSpline(x,v)(np.clip(y,-width,width))
  v=weight(j,x)*np.sum(vals*weights,axis=1)
 return np.sum(CubicSpline(x,v)(np.clip(nodes,-width,width))*weights)

if __name__=='__main__':
 print(sys.version.split()[0],np.__version__,scipy.__version__)
 for U,T in [(0.3,2.),(1.,2.),(3.,1.2)]:
  for N in [8,80,640,2560]:
   for h in ([.1,.01,.001] if N==80 else [.01]):
    st=time.perf_counter();p,b1,b2,M=approximations(N,T,U,h);t1=time.perf_counter()-st
    st=time.perf_counter();r=transfer(N,T,U,h);t2=time.perf_counter()-st
    rr=transfer(N,T,U,h,601,48,12);p2,b3=second_correction(N,T,U,h)
    print(f'U={U:g} T={T:g} N={N:4d} h={h:.3f} Z={rr.real:.12f}{rr.imag:+.12f}i b1={b1:.4g} err1={abs(rr-1):.4g} b2={b2:.4g} err2={abs(rr-p):.4g} b3={b3:.4g} err3={abs(rr-p2):.4g} refdiff={abs(rr-r):.3g} ms={t1*1e3:.3f}/{t2*1e3:.3f}',flush=True)
