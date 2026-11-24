"""Nova 037. Independent floating-point attack, not a global certificate.
Run OPENBLAS_NUM_THREADS=1 python clipped_audit.py.
"""
import numpy as np
from scipy.linalg import eigh, sqrtm, solve, inv
from scipy.optimize import minimize, differential_evolution
from threadpoolctl import threadpool_limits
threadpool_limits(1)

def setup(N):
    d=2/N; w=np.ones(N+1);w[[0,-1]]=.5
    K=(np.diag(np.r_[1,2*np.ones(N-1),1])-np.diag(np.ones(N),1)-np.diag(np.ones(N),-1))/d
    return d,w,K

def saddle(N):
    d,w,K=setup(N);x=np.zeros(N+1);x[0]=-2;p=2+d/2*3*x[0]*np.exp(-x[0]**2/2)
    for j in range(N):
        x[j+1]=x[j]+d*p;p+=d*w[j+1]*3*x[j+1]*np.exp(-x[j+1]**2/2)
    z=x.astype(complex)
    def grad(z,qb,pb):
        f=K@z+3*d*w*z*np.exp(-z*z/2)
        f[0]+=1j*(z[0]+2)+2;f[-1]+=1j*(z[-1]-qb)-pb
        return f
    def hess(z):
        H=K.astype(complex)+np.diag(3*d*w*(1-z*z)*np.exp(-z*z/2));H[0,0]+=1j;H[-1,-1]+=1j
        return H
    for t in np.linspace(0,1,41):
        qb=x[-1]-.6*t;pb=p+.3*t
        for j in range(20):
            step=solve(hess(z),grad(z,qb,pb));z-=step
            if np.linalg.norm(step,np.inf)<1e-13:break
    return z,qb,pb,hess(z),np.linalg.norm(grad(z,qb,pb),np.inf)

def local(N):
    z,qb,pb,H,res=saddle(N);d,w,K=setup(N);s=np.arange(N+1)/N;t=.02+1.4*s*(1-s)
    R=inv(np.asarray(sqrtm(-1j*H),dtype=np.complex128));S=R.imag@inv(R.real)
    C=np.eye(N+1)+1j*S;D=np.diag(1+1j*t);mass=np.diag(d*w)
    return z,qb,pb,H,res,t,S,eigh((C.T@H@C).imag,mass,eigvals_only=True)[0],eigh((D.T@H@D).imag,mass,eigvals_only=True)[0],np.linalg.norm(R.T@(-1j*H)@R-np.eye(N+1),np.inf)

def attack():
    N=32;z,qb,pb,H,res,t,S,*_=local(N);d,w,K=setup(N);G=S-np.diag(t);B=.1
    def phi(Z):
        return .5j*(Z[0]+2)**2+2*Z[0]+np.sum(np.diff(Z)**2)/(2*d)-3*d*np.sum(w*np.exp(-Z*Z/2))+.5j*(Z[-1]-qb)**2-pb*Z[-1]
    base=phi(z).imag
    def fun(v):
        Z=z+v+1j*(t*v+np.clip(G@v,-B,B));return phi(Z).imag-base
    def jac(v):
        gv=G@v;Z=z+v+1j*(t*v+np.clip(gv,-B,B));gr=K@Z+3*d*w*Z*np.exp(-Z*Z/2);gr[0]+=1j*(Z[0]+2)+2;gr[-1]+=1j*(Z[-1]-qb)-pb
        return gr.imag+t*gr.real+G.T@((np.abs(gv)<B)*gr.real)
    T=np.diag(t);D=np.diag(np.r_[1,np.zeros(N-1),1]);A=K@T+T@K+D@(np.eye(N+1)-T*T)
    vals,vec=eigh(A,np.diag(d*w));print('background eig',vals[:5],flush=True)
    rng=np.random.default_rng(1717);best=(0.,np.zeros(N+1)); count=0
    for k in [3,6,10]:
        U=vec[:,:k];bounds=[(-4,4)]*k
        result=differential_evolution(lambda v:fun(U@v),bounds,seed=1717+k,popsize=12,maxiter=250,tol=1e-9,polish=False)
        print('DE',k,result.fun,result.nfev,flush=True)
        rr=minimize(fun,U@result.x,jac=jac,method='BFGS',options={'maxiter':500,'gtol':1e-9});count+=1
        if rr.fun<best[0]:best=(rr.fun,rr.x)
        for j in range(40):
            v=U@(rng.normal(size=k)*rng.choice([.08,.3,1.,3.]))
            rr=minimize(fun,v,jac=jac,method='BFGS',options={'maxiter':300,'gtol':1e-8});count+=1
            if rr.fun<best[0]:best=(rr.fun,rr.x)
        print('best',best[0], 'norm',np.linalg.norm(best[1]),'count',count,flush=True)
    print('bestvector',best[1].tolist(),flush=True)

def diagnostics():
    N=32;z,qb,pb,H,res,t,S,mu,*_=local(N);d,w,K=setup(N);D=np.diag(np.r_[1,np.zeros(N-1),1]);T=np.diag(t);C=np.eye(N+1)+1j*T;B0=(C@(K+1j*D)@C).imag;B=.1;n=N+1
    M=np.diag(d*w);isq=np.diag(1/np.sqrt(d*w));sq=np.diag(np.sqrt(d*w));lam=eigh(B0,M,eigvals_only=True)[0];a=K@z;a[0]+=1j*(z[0]+2)+2;a[-1]+=1j*(z[-1]-qb)-pb
    L=np.linalg.norm(isq@(C.T@a).imag)+B*np.sqrt(n)*np.linalg.norm((K-D@T)@isq,2)
    cc=B*np.sqrt(n)*np.linalg.norm(a.real)+.5*B*B*np.trace(D)+3*d*np.sum(w*(np.exp((abs(z.imag-t*z.real)+B)**2/(2*(1-t*t)))+abs(np.exp(-z*z/2))))
    rad=(L+np.sqrt(L*L+2*lam*cc))/lam
    Ct=np.eye(n)+1j*S;Kinf=max(np.linalg.norm(Ct@isq,axis=1));K2=np.linalg.norm(sq@Ct@isq,2);Ginf=max(np.linalg.norm((S-T)@isq,axis=1))
    M3=max(3*np.exp((abs(zi.imag)+.1)**2/2)*(3*(abs(zi)+.1)+(abs(zi)+.1)**3) for zi in z)
    r=min(B/(2*Ginf),.1/Kinf,1.5*mu/(M3*Kinf*K2*K2))
    print('lam',lam,'L',L,'C',cc,'outer',rad,'local radius',r,'M3',M3,'Kinf',Kinf,'K2',K2,'Ginf',Ginf)

if __name__=='__main__':
    for N in [8,16,32,64,128]:
        o=local(N);print(N, o[4],o[-3:],flush=True)
    attack()
    diagnostics()
