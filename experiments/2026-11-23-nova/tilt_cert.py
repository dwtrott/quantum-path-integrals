import mpmath as mp
from mpmath import iv
import numpy as np
iv.dps=40
A=[[iv.mpc('1.5','.5'),iv.mpc(-1),iv.mpc(0)],[iv.mpc(-1),iv.mpc(2),iv.mpc(-1)],[iv.mpc(0),iv.mpc(-1),iv.mpc('1.5','.5')]]
l=[iv.mpc('1.5','.5'),iv.mpc(0),iv.mpc(0)]
z0=[iv.mpc('-0.3747391002671501','-0.0843695764745657'),iv.mpc('-0.0786082101836452','-0.0173467851927615'),iv.mpc('-0.017680374943359','-0.0018909438209335')]
r=iv.mpf('1e-8'); delta=iv.mpf([-1,1])*r
Z=[z+iv.mpc(delta,delta) for z in z0]
W=lambda z:iv.exp(-z*z/2)
grad=lambda z:[sum(A[i][j]*z[j] for j in range(3))+l[i]+3*z[i]*W(z[i]) for i in range(3)]
hess=lambda z:[[A[i][j]+(3*(1-z[i]*z[i])*W(z[i]) if i==j else 0) for j in range(3)] for i in range(3)]
# Fixed approximate inverse, rounded decimals interpreted as rational intervals.
H0=np.array([[complex(float(v.real.mid),float(v.imag.mid)) for v in row] for row in hess(z0)])
Y=[[iv.mpc(format(v.real,'.12f'),format(v.imag,'.12f')) for v in row] for row in np.linalg.inv(H0)]
F=grad(z0); H=hess(Z)
res=max(abs(sum(Y[i][j]*F[j] for j in range(3))).b for i in range(3))
q=max(sum(abs((1 if i==j else 0)-sum(Y[i][k]*H[k][j] for k in range(3))).b for j in range(3)).b for i in range(3))
assert q<1 and res+q*r<r
print('root residual upper',res,'contraction upper',q)
limits=['-.34','-.06','-.015']; tilts=['0','.05','0']
for z,t0,L in zip(Z,tilts,limits):
 t=iv.mpf(t0); c=iv.mpc(1,t); L=iv.mpf(L); R=64
 k=z.imag-t*z.real; M=iv.exp(k*k/(2*(1-t*t)))
 tail=-6*((M+abs(W(z)))/R**2+abs(z*W(z)*c)/R)
 assert tail.a>L, (tail,L)
 near=iv.mpf(['-.0625','.0625']); v=z+c*near
 nearbound=(3*(1-v*v)*W(v)*c*c).imag
 assert nearbound.a>L,(nearbound,L)
 stack=[(-R,-.0625),(.0625,R)]; count=0; depth=0
 while stack:
    a,b=stack.pop(); s=iv.mpf([a,b]); v=z+c*s
    val=(-6*(W(v)-W(z)+z*W(z)*c*s).imag/(s*s))
    if val.a>L: count+=1;continue
    if b-a<1e-6: raise RuntimeError((a,b,val,L))
    mid=(a+b)/2;stack.extend([(a,mid),(mid,b)])
 print('scalar bound',L,'panels',count,'near',nearbound,'tail',tail)
# B+diag(L) at t=(0,.05,0): diagonal (.16,.14,.485), offdiag -.05.
# Gershgorin >= .04 I.
