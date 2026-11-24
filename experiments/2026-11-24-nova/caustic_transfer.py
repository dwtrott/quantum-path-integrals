import math,time,sys
import numpy as np
import mpmath as mp
from mpmath import iv
iv.dps=40
I=iv.mpf

def up(x): return math.nextafter(float(x.b),math.inf)
def mid(x): return float(x.mid)
def theta_mass(a):
    v=iv.exp(-a); return 1+2*v/(1-v)
def setup(N,eta,M,dxs):
    h=I('1e-5'); d=iv.sqrt(2)*iv.sin((iv.pi+I(eta))/(2*N)); T=N*d; n=N-1
    ang=iv.pi/5; s=iv.sin(2*ang); c=iv.cos(2*ang); k=I('.004'); K=12; dx=I(dxs); r=I('.0017')
    B=[2/d*(iv.cos((iv.pi+I(eta))/N)-iv.cos(j*iv.pi/N)) for j in range(1,N)]
    P=[s*x+k*d for x in B]; assert min(x.a for x in P)>0
    A0=(c*(1-I(1)/K**2)-2*s/K)/2; D0=(s*(1-I(1)/K**2)-2*c/K)/2
    assert A0.a>0 and D0.a>0
    eps=k/(2*D0); C0=(D0/A0)*(1-iv.sqrt(1-eps))**2
    Cshift=k*K*K*r*r/2+(K*K+1)**2*r**4*iv.exp((K*K+1)*r*r/2)/8
    Z=iv.sqrt(4*iv.pi*h/N)*iv.exp((sum(iv.ln(x) for x in B[1:])-sum(iv.ln(x) for x in P))/2)/iv.sqrt(iv.sqrt(h))
    mass=Z*iv.exp(C0*T/h)
    L=max((s*x+c*c*x*x/y).b for x,y in zip(B,P))
    v=max(sum(2*h*iv.sin(j*l*iv.pi/N)**2/(N*p) for l,p in enumerate(P,1)).b for j in range(1,N))
    a=2*iv.pi*r/dx-L*r*r/(2*h)-d/h*Cshift
    assert a.a>0
    alias=mass*(theta_mass(a)**n-1)
    aP=2*iv.pi**2*h/(dx**2*max(x.b for x in P)); box=mass*theta_mass(aP)**n*2*n*iv.exp(-(M*dx)**2/(2*v))
    b=s/d; diag=2*b-2*s*d+k*d; prec=[diag]
    for j in range(n-1): prec.append(diag-b*b/prec[-1])
    prec=prec[::-1]; sig2=[h/p for p in prec]; assert min(p.a for p in prec)>0
    th=theta_mass(2*iv.pi**2*min(x.a for x in sig2)/dx**2)
    band=mass*th**n*2*n*iv.exp(-50)
    print('N',N,'eta',eta,'M',M,'dx',dxs,'mass',up(mass),'alias',up(alias),'box',up(box),'band',up(band),'theta',up(th),flush=True)
    return locals()

def run(N,eta,M,dxs,validated=False):
    t0=time.perf_counter(); p=setup(N,eta,M,dxs); globals0=p
    h,d,s,c,k,dx,Z,mass,C0,prec,sig2,b,th,ang=[p[x] for x in ('h','d','s','c','k','dx','Z','mass','C0','prec','sig2','b','th','ang')]; n=N-1
    idx=np.arange(-M,M+1,dtype=np.int64); D=mid(dx); yy=idx*D
    def c_mid(z): return complex(mid(z.real),mid(z.imag))
    def c_err(z,zz): return up(abs(z-iv.mpc(zz.real,zz.imag)))
    def site(j):
        x=j*dx; q=x*x/2; w=iv.mpc(0,d/h)*(iv.exp(-iv.mpc(c,s)*q)-1+iv.mpc(c,s)*q)
        return iv.exp(w+k*d*q/h-iv.mpc(0,2*c*d*q/h))
    def edge(j): return iv.exp(iv.mpc(0,c*(j*dx)**2/(2*h*d)))
    R0=[site(j) for j in range(M+1)]; E0=[edge(j) for j in range(2*M+1)]
    rv=np.array([c_mid(x) for x in R0]); ev=np.array([c_mid(x) for x in E0]); er=max(c_err(z,v0) for z,v0 in zip(R0,rv)); ee=max(c_err(z,v0) for z,v0 in zip(E0,ev)); print('factor errors',er,ee,flush=True)
    F=ev[np.abs(idx)]; fac=rv[np.abs(idx)]; u=2**-53; rowerr=I(0); rowmax=iv.exp(C0*d/h)*th
    E=I(ee); Bnd=I(1); maxW=0
    for j in range(n-1,-1,-1):
        delta=1-b/prec[j]; ac=dx**2/(2*sig2[j]); pc=dx/iv.sqrt(2*iv.pi*sig2[j]); dl,aa,pp=map(mid,(delta,ac,pc));
        ii=idx if j else np.array([0],dtype=np.int64)
        centers=np.rint((1-dl)*ii).astype(np.int64); W=math.ceil(10*math.sqrt(up(sig2[j]))/float(dx.a))+3; maxW=max(maxW,W)
        kk=centers[:,None]+np.arange(-W,W+1); dif=kk-ii[:,None]; good=(kk>=-M)&(kk<=M)
        assert (I(W)-I(2))*dx.a > 10*iv.sqrt(sig2[j]).b
        assert up(M*abs(I(dl)-delta)+I(4)*I(u)*M*(1+abs(I(dl)))) < .01
        DD=dif+dl*ii[:,None]; tt=aa*DD*DD
        de=I(M)*abs(I(dl)-delta)+I(2)*I(u)*abs(I(dl))*M+I(2)*I(u)*(I(int(np.max(np.abs(dif))))+abs(I(dl))*M)
        Dmax=I(float(np.max(np.abs(DD))))+de
        te=abs(I(aa)-ac)*Dmax**2+abs(I(aa))*(2*Dmax*de+de**2)+I(4)*I(u)*abs(I(aa))*Dmax**2
        assert (I(float(tt.max()))+te).b < 128 and tt.min()>=0
        if validated:
            x=tt/128; w=np.full_like(x,1/math.factorial(16))
            for l in range(15,-1,-1): w=1/math.factorial(l)-x*w
            for _ in range(7): w=w*w
        else: w=np.exp(-tt)
        weight=pp*w; jj=np.clip(kk+M,0,2*M); val=F[jj]*fac[jj]*ev[np.minimum(np.abs(dif),2*M)]; val*=weight; val[~good]=0
        F=np.sum(val,axis=1)
        ce=abs(I(pp)-pc)+abs(I(pp))*(te+I('4e-12'))+I(2)*I(u)*abs(I(pp));
        Krow=2*W+1; Rmax=iv.exp(C0*d/h); rowerr=Krow*ce*(Rmax+I(er))*(1+I(ee)) + th*(I(er)+Rmax*I(ee)+I(er)*I(ee))
        rowerr+=I(32*Krow)*I(u)*(rowmax+rowerr)+I('1e-280')
        E=(rowmax+rowerr)*E+rowerr*Bnd; Bnd*=rowmax
    pref=Z*iv.exp(iv.mpc(0,(ang-iv.pi/4)*n)); pf=c_mid(pref); pe=c_err(pref,pf)
    out=pf*F[0]; roundbound=abs(pref)*E+I(pe)*(Bnd+E)+I(16)*I(u)*abs(pref)*(Bnd+E)
    print('value',out,'round_bound',up(roundbound) if validated else 'not validated','maxW',maxW,'sec',time.perf_counter()-t0,flush=True)
    print('total',up(p['alias']+p['box']+p['band']+roundbound) if validated else 'not validated',flush=True)
    return out
if __name__=='__main__':
    print(sys.version,np.__version__,mp.__version__,flush=True)
    for N,eta,M,dx in [(4,'0',1600,'.00022'),(80,'-.001',1800,'.00022'),(80,'.001',1800,'.00022')]: run(N,eta,M,dx,True)
