import numpy as np
from scipy.linalg import eigh_tridiagonal
from mpmath import iv
import mpmath, time
iv.dps=35
I=iv.mpf; C=iv.mpc
W=lambda z:iv.exp(-z*z/2)
# interval endpoints extracted without conversion to binary64
# choose lower ell using float guide only; all accepted tests interval comparisons

def certify(N):
 d=I(2)/N;t=I(1)/10;c=C(1,t);x=I(-2);p=I(2)+d/2*3*x*W(x);xs=[x]
 for j in range(N):
  x=x+d*p;xs.append(x);p=p+d*(I(1)/2 if j==N-1 else 1)*3*x*W(x)
 weights=np.ones(N+1);weights[[0,-1]]=.5
 us=np.r_[-np.geomspace(128,1e-4,6000),np.geomspace(1e-4,128,6000)]
 ell=[];counts=[]
 for X in xs:
  z=float(X.mid);w=np.exp(-z*z/2);rz=2*np.imag(-3*(np.exp(-(z+(1+.1j)*us)**2/2)-w+(1+.1j)*us*z*w))/us**2
  lower=I(int(np.floor((rz.min()-.035)*1000)))/1000
  # near zero Taylor Hessian enclosure
  r=1/128;R=128
  def interval(a,b):return I([a,b])
  Z=C(X,0);w=W(Z);wp=-Z*w
  def hess(a,b):
   v=Z+c*interval(a,b);return (-3*c*c*(v*v-1)*W(v)).imag
  assert hess(-r,r)>lower,(N,z,'near',hess(-r,r),lower)
  M=iv.exp((t*X)**2/(2*(1-t*t)))
  tail=-6*((M+abs(w))/R**2+abs(c*wp)/R)
  assert tail>lower,(N,z,'tail',tail,lower)
  stack=[(-R,-r),(r,R)];n=0
  while stack:
   a,b=stack.pop();u=interval(a,b)
   val=2*(-3*(W(Z+c*u)-w-c*u*wp)).imag/(u*u)
   if val>lower:n+=1;continue
   assert b-a>2**-32,(N,z,a,b,val,lower)
   m=(a+b)/2;stack.extend([(a,m),(m,b)])
  ell.append(lower);counts.append(n)
 # Sylvester via interval LDL on shifted pencil.  mu=.12 for both
 mu=I(12)/100;off=-2*t/d;P=None
 for j,e in enumerate(ell):
  end=j in (0,N);a=(2*t/d)*(1 if end else 2)+(1-t*t if end else 0)+d*(I(1)/2 if end else 1)*(e-mu)
  P=a if P is None else a-off*off/P
  assert P>0,(N,j,P)
 # actual eig lower approx diagnostic
 dd=2/N;wf=weights;diag=2*.1/dd*np.r_[1,np.full(N-1,2),1]+.99*np.r_[1,np.zeros(N-1),1]+dd*wf*np.array([float(e.mid) for e in ell]);offf=-2*.1/dd/np.sqrt((dd*wf[:-1])*(dd*wf[1:]));ev=eigh_tridiagonal(diag/(dd*wf),offf,select='i',select_range=(0,0))[0][0]
 print(N,'panels',sum(counts),'min ell',min(float(e.mid) for e in ell),'eig diagnostic',ev,'last pivot',P,'x width',xs[-1].delta,flush=True)
if __name__=='__main__':
 for N in [16,64]:certify(N)
