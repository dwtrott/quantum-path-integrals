"""Independent floating-point diagnostics for Sparks 038. Not interval verification.
OPENBLAS_NUM_THREADS=1 python structured_audit.py
"""
import numpy as np
import scipy
from scipy.linalg import eigh, solve
from scipy.optimize import differential_evolution, minimize
from threadpoolctl import threadpool_limits
threadpool_limits(1)

def setup(N):
    d=2/N; m=d*np.r_[.5,np.ones(N-1),.5]
    E=np.diff(np.eye(N+1),axis=0); K=E.T@E/d
    D=np.diag(np.r_[1,np.zeros(N-1),1]);t=.02+1.4*(np.arange(N+1)/N)*(1-np.arange(N+1)/N)
    x=np.zeros(N+1);x[0]=-2;p=2+1.5*d*x[0]*np.exp(-x[0]**2/2)
    for j in range(N):
        x[j+1]=x[j]+d*p;p+=3*m[j+1]*x[j+1]*np.exp(-x[j+1]**2/2)
    z=x.astype(complex)
    def grad(Z,qb,pb):
        G=K@Z+3*m*Z*np.exp(-Z*Z/2);G[0]+=1j*(Z[0]+2)+2;G[-1]+=1j*(Z[-1]-qb)-pb;return G
    for u in np.linspace(0,1,41):
        qb=x[-1]-.6*u;pb=p+.3*u
        for j in range(30):
            H=K+1j*D+np.diag(3*m*(1-z*z)*np.exp(-z*z/2));step=solve(H,grad(z,qb,pb));z-=step
            if np.max(abs(step))<1e-13:break
    def phi(Z):return .5j*(Z[0]+2)**2+2*Z[0]+np.sum(np.diff(Z)**2)/(2*d)-3*np.sum(m*np.exp(-Z*Z/2))+.5j*(Z[-1]-qb)**2-pb*Z[-1]
    return d,m,E,K,D,t,z,phi,lambda Z:grad(Z,qb,pb)

def run(N,attack=False):
    d,m,E,K,D,t,z,phi,grad=setup(N); M=np.diag(m);T=np.diag(t);sigma=.8291;a=-1.6247;B=.4;V=np.ones(N+1)/np.sqrt(2);P=V*m
    A=K@T+T@K+D@(np.eye(N+1)-T*T);vals,U=eigh(A,M)
    H=K+1j*D+np.diag(3*m*(1-z*z)*np.exp(-z*z/2));C=np.eye(N+1)+1j*(sigma*np.eye(N+1)+a*np.outer(V,P));mu=eigh((C.T@H@C).imag,M,eigvals_only=True)[0]
    base=phi(z).imag
    def mapv(v,b=B):
        u=(sigma-t)*v;q=P@v;return z+v+1j*(t*v+np.clip(u,-b,b)+a*V*np.clip(q,-b,b)),u,q
    def fun(v,b=B):return phi(mapv(v,b)[0]).imag-base
    def jac(v,b=B):
        Z,u,q=mapv(v,b);g=grad(Z);return g.imag+(t+(abs(u)<b)*(sigma-t))*g.real+a*(abs(q)<b)*P*(V@g.real)
    # Exterior: reserve eta of average-tilt kinetic energy to absorb clip variation.
    # At each edge |clip((sigma-t_j)x)-clip((sigma-t_i)x)| <= B |dt|/(sigma-max(t_i,t_j)).
    av=(t[1:]+t[:-1])/2; R=E.T@np.diag(av/d)@E
    c=B*abs(np.diff(t))/(sigma-np.maximum(t[1:],t[:-1]));beta=B*(1+abs(a)/np.sqrt(2));gq=-3*m*z*np.exp(-z*z/2)
    # endpoint linear from -t v (b+c), absorb in mass-independent norm via Riesz with B_remaining.
    ell=(gq.imag+t*gq.real).copy() # actually Im((1+it)gq) = imag + t real
    ell_abs=abs(ell);ell_abs[[0,-1]]+=beta*t[[0,-1]]
    E0=beta*np.sum(abs(gq.real))+beta**2+3*np.sum(m*(np.exp((abs(z.imag-t*z.real)+beta)**2/(2*(1-t*t)))+abs(np.exp(-z*z/2))))
    # sharpen endpoint/linear bound using A_remaining^{-1}, sign-independent bound by abs inverse
    best=None
    for eta in [.1,.2,.3,.4,.5,.6,.7,.8,.9]:
        Ar=A-2*eta*R;lam=eigh(Ar,M,eigvals_only=True)[0]
        if lam<=0:continue
        E1=E0+np.sum(c*c/(4*eta*av*d));L=np.sqrt(ell_abs@abs(np.linalg.inv(Ar))@ell_abs)
        # h >= .5*||v||_Ar^2 - L*||v||_Ar - E1
        rad=(L+np.sqrt(L*L+2*E1))/np.sqrt(lam)
        if best is None or rad<best[0]:best=(rad,eta,lam,L,E1)
    print('N',N,'versions',np.__version__,scipy.__version__,'res',max(abs(grad(z))),'local',mu,'background',vals[0],'exterior',best,flush=True)
    if not attack:return
    rng=np.random.default_rng(2206);best=(0.,np.zeros(N+1));nstart=0
    for k in [3,6,10]:
        Vsub=U[:,:k];r=differential_evolution(lambda w:fun(Vsub@w),[(-4,4)]*k,seed=2206+k,popsize=12,maxiter=250,tol=1e-9,polish=False)
        print('DE',k,r.fun,r.nfev,flush=True)
        for w in [r.x]+[rng.normal(size=k)*rng.choice([.08,.3,1.,3.]) for _ in range(40)]:
            rr=minimize(fun,Vsub@w,jac=jac,method='BFGS',options={'maxiter':400,'gtol':1e-8});nstart+=1
            if rr.fun<best[0]:best=(rr.fun,rr.x)
        print('BFGS',nstart,best[0],np.linalg.norm(best[1]),flush=True)
    # positive control: narrower clip known to fail. Independent search on same modes.
    r=differential_evolution(lambda w:fun(U[:,:6]@w,.1),[(-4,4)]*6,seed=2207,popsize=12,maxiter=250,polish=False)
    rr=minimize(lambda v:fun(v,.1),U[:,:6]@r.x,jac=lambda v:jac(v,.1),method='BFGS',options={'maxiter':400,'gtol':1e-8})
    print('negative control B=.1',rr.fun,'v',rr.x.tolist(),flush=True)

if __name__=='__main__':
    for N in [16,32,64,128,256]:run(N,N==32)
