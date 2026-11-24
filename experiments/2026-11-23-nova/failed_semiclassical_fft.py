# Nova, 2026-11-23. NEGATIVE RESULT, not a usable integrator.
# Run: OPENBLAS_NUM_THREADS=1 python failed_semiclassical_fft.py
# Multiplication by the positive quadratic weights amplifies FFT noise.
import numpy as np
from scipy.signal import fftconvolve

def green_diag(N,T,U,m=1.,omega=1.):
 d=T/N;nu=np.sqrt(omega**2+U/m);theta=2*np.arcsin(d*nu/2)
 assert N*theta<np.pi
 j=np.arange(1,N);return d*np.sin(j*theta)*np.sin((N-j)*theta)/(m*np.sin(theta)*np.sin(N*theta))

def transfer(N,T,U,h,points_per_step=12,width=10):
 d=T/N;M=green_diag(N,T,U).max();dx=np.sqrt(d)/points_per_step;K=int(np.ceil(width*np.sqrt(M)/dx));x=np.arange(-K,K+1)*dx
 offsets=np.arange(-2*K,2*K+1)*dx;kernel=np.exp(-offsets**2/(2*d))*dx/np.sqrt(2*np.pi*d)
 a=np.exp(-x*x/(2*d))/np.sqrt(2*np.pi*d);v=h*x*x/2
 f=np.exp(1j*U*d/h*(np.expm1(-1j*v)+1j*v));w=np.exp(d*(1+U)*x*x/2)
 arr=np.stack([a.astype(complex),a.astype(complex)])
 for _ in range(N-1):arr=fftconvolve(arr*np.stack([w*f,w]),kernel[None,:],mode='same',axes=1)
 return arr[0,K]/arr[1,K],2*K+1

if __name__=='__main__':
 for N in (8,80):
  print(N,transfer(N,2,.3,.01),transfer(N,2,.3,.01,16,11))
