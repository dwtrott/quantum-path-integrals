import numpy as np, math
from scipy.optimize import differential_evolution

def search(N,eta,h=1e-5,target=.01):
 n=N-1; d=math.sqrt(2)*math.sin((math.pi+eta)/(2*N)); T=N*d
 B=2/d*(math.cos((math.pi+eta)/N)-np.cos(np.arange(1,N)*math.pi/N))
 def calc(x,ret=False):
  ang,k,K,r,dx=x; s=math.sin(2*ang); c=math.cos(2*ang); P=s*B+k*d
  A=(c*(1-1/K**2)-2*s/K)/2; D=(s*(1-1/K**2)-2*c/K)/2
  if min(P)<=0 or A<=0 or D<=k/2:return 1e30 if not ret else None
  def C(A,D):return (D/A)*(k/(2*D)/(1+math.sqrt(1-k/(2*D))))**2
  C0=C(c/2,s/2); Cs=max(C(A,D),k*K*K*r*r/2+(K*K+1)**2*r**4*math.exp(min(100,(K*K+1)*r*r/2))/8)-C0; Cs=max(0,Cs)
  logmass=.5*(math.log(4*math.pi*h/N)-.5*math.log(h)+np.log(B[1:]).sum()-np.log(P).sum())+C0*T/h
  L=max(s*B+c*c*B*B/P); a=2*math.pi*r/dx-L*r*r/(2*h)-d/h*Cs
  if a<=0 or logmass>100:return 1e30 if not ret else None
  # stable log Theta = log1p(2/(exp(a)-1))
  def lt(a): return math.log1p(2/math.expm1(a)) if a<700 else 0
  alias=math.exp(logmass)*math.expm1(n*lt(a)) if n*lt(a)<100 else 1e30
  b=s/d; diag=2*b-2*s*d+k*d; prec=[diag]
  for j in range(n-1):prec.append(diag-b*b/prec[-1])
  prec=np.array(prec[::-1]); sig=h/prec; v=0; vmax=0
  for p,ss in zip(prec,sig):v=ss+(b/p)**2*v;vmax=max(vmax,v)
  lP=lt(2*math.pi**2*h/(dx*dx*max(P))); lth=lt(2*math.pi**2*min(sig)/dx**2)
  M=math.ceil(math.sqrt(2*vmax*(logmass+n*lP+math.log(2*n/(target*.3))))/dx)
  W=math.ceil(10*math.sqrt(max(sig))/dx)+3
  # approximate conservative budget of 4e-12 absolute exponential error
  if logmass+n*lth>100:return 1e30 if not ret else None
  roundest=math.exp(logmass+n*lth)*n*(2*W+1)*max(dx/np.sqrt(2*math.pi*sig))*4.1e-12
  band=math.exp(logmass+n*lth)*2*n*math.exp(-50)
  cost=(2*M+1)*(2*W+1)*n
  penalty=max(0,math.log(max(1e-300,alias+roundest+band)/(target*.7)))
  if ret:return dict(ang=ang,k=k,K=K,r=r,dx=dx,M=M,W=W,alias=alias,roundest=roundest,band=band,mass=math.exp(logmass),cost=cost)
  return math.log(cost)+10*penalty
 bounds=[(.54,.784),(.00201,.014),(4,150),(.0001,.006),(.000025,.001)]
 res=differential_evolution(calc,bounds,seed=17,popsize=10,maxiter=110,tol=.001); print(N,eta,calc(res.x,True),flush=True)
if __name__=='__main__':
 for N in [80,320,1280]:search(N,.001)
