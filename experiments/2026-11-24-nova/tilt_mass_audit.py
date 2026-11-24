"""Independent shooting and Hessian mass-normalization audit; floating point only."""
import numpy as np
import scipy
from scipy.optimize import root, minimize
from scipy.linalg import eigh_tridiagonal
W=lambda z:np.exp(-z*z/2)
def shoot(z,N):
 d=2/N; p=2+1j*(z+2)+d/2*3*z*W(z); xs=[z]
 for j in range(N):
  z=z+d*p; xs.append(z); p=p+d*(.5 if j==N-1 else 1)*3*z*W(z)
 return np.array(xs),p

def saddle(N,dq,dp):
 x,p=shoot(-2,N); qb=x[-1].real; pb=p.real; z=-2+0j
 for s in np.linspace(0,1,21)[1:]:
  def f(a):
   y,v=shoot(complex(*a),N); r=v-(pb+s*dp)+1j*(y[-1]-qb-s*dq); return [r.real,r.imag]
  rr=root(f,[z.real,z.imag],tol=1e-11); assert np.linalg.norm(rr.fun)<1e-10; z=complex(*rr.x)
 return shoot(z,N)[0]

def eigdiag(z,t,mass):
 N=len(z)-1;d=2/N;w=np.ones(N+1);w[[0,-1]]=.5;D=np.r_[1,np.zeros(N-1),1];K=np.r_[1,np.full(N-1,2),1]/d
 c=1+1j*t;diag=2*t*K+(1-t*t)*D+np.imag(-3*d*w*c*c*(z*z-1)*W(z));off=-(t[:-1]+t[1:])/d
 if mass:diag=diag/(d*w);off=off/(d*np.sqrt(w[:-1]*w[1:]))
 return eigh_tridiagonal(diag,off,select='i',select_range=(0,0))[0][0]

def main():
 print(np.__version__,scipy.__version__,flush=True)
 rng=np.random.default_rng(2)
 print('N MorseIndex realH_min raw_at_mass_opt mass_at_mass_opt raw_at_raw_opt mass_at_raw_opt constant_t_classical_raw constant_t_classical_mass',flush=True)
 for N in [8,16,32,64]:
  z=saddle(N,-.6,.3);d=2/N;w=np.ones(N+1);w[[0,-1]]=.5;diag=np.r_[1,np.full(N-1,2),1]/d-3*d*w*np.real((z*z-1)*W(z));ev=eigh_tridiagonal(diag,np.full(N,-1/d),eigvals_only=True)
  opts=[]
  for mass in [True,False]:
   best=None
   for start in [np.full(N+1,.2),.02+1.4*np.linspace(0,1,N+1)*(1-np.linspace(0,1,N+1)),rng.uniform(.0,.3,N+1)]:
    rr=minimize(lambda t:-eigdiag(z,t,mass),start,method='L-BFGS-B',bounds=[(-.95,.95)]*(N+1),options={'maxiter':300,'ftol':1e-12,'gtol':1e-8})
    if best is None or rr.fun<best.fun:best=rr
   opts.append(best.x)
  a,b=opts;zc=saddle(N,0,0);tc=np.full(N+1,.1)
  print(N,sum(ev<0),ev[0],eigdiag(z,a,False),eigdiag(z,a,True),eigdiag(z,b,False),eigdiag(z,b,True),eigdiag(zc,tc,False),eigdiag(zc,tc,True),flush=True)
if __name__=='__main__':main()
