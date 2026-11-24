import numpy as np
from scipy.optimize import root, minimize_scalar
from scipy.linalg import eigvals
A=np.array([[1.5+.5j,-1,0],[-1,2,-1],[0,-1,1.5+.5j]])
l=np.array([1.5+.5j,0,0]);g=3
W=lambda z:np.exp(-z*z/2)
grad=lambda z:A@z+l+g*z*W(z)
H=lambda z:A+np.diag(g*(1-z*z)*W(z))
phi=lambda z:.5*z@A@z+l@z+2+(-1+1j)/4-g*sum(W(z))
rng=np.random.default_rng(42); roots=[]
with np.errstate(all='ignore'):
 for s in range(1500):
    r=root(lambda x:np.r_[grad(x[:3]+1j*x[3:]).real,grad(x[:3]+1j*x[3:]).imag],np.r_[rng.uniform(-4,4,3),rng.uniform(-1.5,1.5,3)])
    z=r.x[:3]+1j*r.x[3:]
    if np.linalg.norm(grad(z))<1e-8 and np.max(abs(z.imag))<3 and not any(np.linalg.norm(z-w)<1e-6 for w in roots):roots.append(z)
def reference(h,n=600,m=1600):
    x=np.linspace(-1-10*np.sqrt(h),-1+10*np.sqrt(h),n+1); y=np.linspace(-10*np.sqrt(h),10*np.sqrt(h),n+1)
    u=np.linspace(-8,8,m+1); z=(1+.3j)*u
    w1=np.ones(n+1);w1[[0,-1]]=.5; w2=np.ones(m+1);w2[[0,-1]]=.5
    left=np.exp(1j*((x+2)**2/2+(-1+1j)*(x+1)**2/4-g*W(x))[None,:]/h)*np.exp(1j*(z[:,None]-x[None,:])**2/(2*h))
    right=np.exp(1j*(y*y/2+(-1+1j)*y*y/4-g*W(y))[None,:]/h)*np.exp(1j*(z[:,None]-y[None,:])**2/(2*h))
    return (1+.3j)*np.sum(w2*(left@w1)*(right@w1)*np.exp(-1j*g*W(z)/h))*(x[1]-x[0])*(y[1]-y[0])*(u[1]-u[0])
I=reference(.05); I1=reference(.1)
def lead(z,h,t=np.array([0,.05,0])):
    C=np.diag(1+1j*t); mat=-1j*C@H(z)@C
    return np.linalg.det(C)*(2*np.pi*h)**1.5/np.prod(np.sqrt(eigvals(mat)))*np.exp(1j*phi(z)/h)
rs=sorted([z for z in roots if phi(z).imag>=0],key=lambda z:phi(z).imag)
best=min(rs,key=lambda z:min(abs(I-lead(z,.05)),abs(I+lead(z,.05))))
print('roots',len(roots),len(rs),'I',I,'I1',I1,'best',best,'action',phi(best),'lowest',rs[0],phi(rs[0]),'err',abs(I-lead(best,.05))/abs(I))
print('reference refinement',reference(.05,1200,3200),reference(.1,1200,3200))
for t in [np.array([0,.05,0]),np.array([0,.1,0]),np.array([.1,.1,.1]),np.array([-.1,.1,.1])]:
 C=np.diag(1+1j*t); B=(C@A@C).imag
 u=np.linspace(-15,15,60001);u=u[abs(u)>1e-4]
 vals=[]
 for z,c in zip(best,np.diag(C)):
    r=-2*g*(W(z+c*u)-W(z)+z*W(z)*c*u).imag/u**2
    vals.append(min(r.min(),(g*(1-z*z)*W(z)*c*c).imag))
 print('tilt',t,'B',B,'remainder minima',vals,'eigs',np.linalg.eigvalsh(B+np.diag(vals)))
np.set_printoptions(precision=16); print('best precise',best)
for t in [np.array([0,.05,0]),np.array([0,.1,0]),np.array([.1,.1,.1])]:
 C=np.diag(1+1j*t); B=(C@A@C).imag
 u=np.linspace(-15,15,60001);vals=[(g*(1-(z+c*u)**2)*W(z+c*u)*c*c).imag.min() for z,c in zip(best,np.diag(C))];print('Hessian',t,vals,np.linalg.eigvalsh(B+np.diag(vals)))
