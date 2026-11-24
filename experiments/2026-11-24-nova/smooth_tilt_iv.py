"""Interval smooth tilt certificate at exact discrete classical target."""
import numpy as np
import mpmath,time
from mpmath import iv
iv.dps=35
I=iv.mpf;C=iv.mpc;W=lambda z:iv.exp(-z*z/2)
def certify(N=64):
 d=I(2)/N;x=I(-2);p=I(2)+d/2*3*x*W(x);xs=[x]
 for j in range(N):
  x=x+d*p;xs.append(x);p=p+d*(I(1)/2 if j==N-1 else 1)*3*x*W(x)
 ts=[I(1)/50+I(7)/5*(I(j)/N)*(1-I(j)/N) for j in range(N+1)]
 us=np.r_[-np.geomspace(128,1e-4,6000),np.geomspace(1e-4,128,6000)];ell=[];count=0
 for j,(X,t) in enumerate(zip(xs,ts)):
  z=float(X.mid);tf=float(t.mid);w0=np.exp(-z*z/2);c0=1+1j*tf
  rz=2*np.imag(-3*(np.exp(-(z+c0*us)**2/2)-w0+c0*us*z*w0))/us**2
  lower=I(int(np.floor((rz.min()-.035)*1000)))/1000
  c=C(1,t);Z=C(X,0);w=W(Z);wp=-Z*w;r=1/128;R=128
  v=Z+c*I([-r,r]);hess=(-3*c*c*(v*v-1)*W(v)).imag
  assert hess.a>lower.b,(j,'near',hess,lower)
  M=iv.exp((t*X)**2/(2*(1-t*t)));tail=-6*((M+abs(w))/R**2+abs(c*wp)/R)
  assert tail.a>lower.b,(j,'tail',tail,lower)
  stack=[(-R,-r),(r,R)]
  while stack:
   a,b=stack.pop();u=I([a,b]);val=2*(-3*(W(Z+c*u)-w-c*u*wp)).imag/(u*u)
   if val.a>lower.b:count+=1;continue
   assert b-a>2**-32,(j,a,b,val,lower)
   mid=(a+b)/2;stack.extend([(a,mid),(mid,b)])
  ell.append(lower)
 mu=I(35)/100;P=None;small=None
 for j,(t,e) in enumerate(zip(ts,ell)):
  end=j in (0,N);a=2*t/d*(1 if end else 2)+(1-t*t if end else 0)+d*(I(1)/2 if end else 1)*(e-mu)
  off=-(ts[j-1]+t)/d if j else None
  P=a if P is None else a-off*off/P
  assert P.a>0,(j,P)
  if small is None or P.a<small:small=P.a
 print('interval',N,'mu',mu,'panels',count,'min ell',min(float(e.mid) for e in ell),'min pivot lower',small,'last pivot',P,flush=True)
if __name__=='__main__':
 print(np.__version__,mpmath.__version__,flush=True);certify()
