"""Improved Gaussian envelope for PR 34's full lattice; no modes removed."""
import time,sys
import caustic_transfer_folded as old
from caustic_transfer_folded import I,iv,up,theta_mass
from full_poly_budget import full as soft_run

def setup(N,eta,M,dxs,ang_s,k_s,K,r_s):
    h=I('1e-5'); d=iv.sqrt(2)*iv.sin((iv.pi+I(eta))/(2*N)); T=N*d; n=N-1
    ang=I(ang_s); s=iv.sin(2*ang); c=iv.cos(2*ang); k=I(k_s); dx=I(dxs); r=I(r_s)
    B=[2/d*(iv.cos((iv.pi+I(eta))/N)-iv.cos(j*iv.pi/N)) for j in range(1,N)]
    P=[s*x+k*d for x in B]; assert min(x.a for x in P)>0
    A0=(c*(1-I(1)/K**2)-2*s/K)/2; D0=(s*(1-I(1)/K**2)-2*c/K)/2
    assert A0.a>0 and D0.a>0
    def C(A,D):
        e=k/(2*D); assert e.a>0 and e.b<1
        return (D/A)*(e/(1+iv.sqrt(1-e)))**2
    C0=C(c/2,s/2)
    Ck=C(A0,D0); Cr=k*K*K*r*r/2+(K*K+1)**2*r**4*iv.exp((K*K+1)*r*r/2)/8
    Cshift=max(I(0), (Ck-C0).b, (Cr-C0).b)
    Z=iv.sqrt(4*iv.pi*h/N)*iv.exp((sum(iv.ln(x) for x in B[1:])-sum(iv.ln(x) for x in P))/2)/iv.sqrt(iv.sqrt(h))
    mass=Z*iv.exp(C0*T/h)
    L=max((s*x+c*c*x*x/y).b for x,y in zip(B,P))
    a=2*iv.pi*r/dx-L*r*r/(2*h)-d/h*Cshift; assert a.a>0
    alias=mass*(theta_mass(a)**n-1)
    b=s/d; diag=2*b-2*s*d+k*d; prec=[diag]
    for j in range(n-1):
        assert prec[-1].a>0
        prec.append(diag-b*b/prec[-1])
    prec=prec[::-1]; assert min(p.a for p in prec)>0
    sig2=[h/p for p in prec]
    var=I(0); v=I(0); vv=[]
    for p,ss in zip(prec,sig2):
        var=ss+(b/p)**2*var; vv.append(var); v=max(v,var.b)
    aP=2*iv.pi**2*h/(dx**2*max(x.b for x in P)); box=mass*theta_mass(aP)**n*2*n*iv.exp(-(M*dx)**2/(2*v))
    th=theta_mass(2*iv.pi**2*min(x.a for x in sig2)/dx**2)
    band=mass*th**n*2*n*iv.exp(-50)
    print('optimized N',N,'eta',eta,'M',M,'dx',dxs,'angle',ang_s,'k',k_s,'K',K,'r',r_s,'mass',up(mass),'alias',up(alias),'box',up(box),'band',up(band),'theta',up(th),flush=True)
    return locals()

CASES={80:(1075,'.0002496','.61447','.0065895',48,'.00162085'),320:(2773,'.00012245','.68714','.0048462',93,'.00086119'),1280:(7725,'.000059475','.73714','.0033275',83,'.00046067')}

def run(N,eta='.001'):
    M,dx,*par=CASES[N]
    old.setup=lambda N,eta,M,dxs:setup(N,eta,M,dxs,*par)
    return old.run(N,eta,M,dx,True)

def check(N):
    M,dx,*par=CASES[N]; p=setup(N,'.001',M,dx,*par)
    for j in [1,N//2,N-1]:
        spec=sum(2*p['h']*iv.sin(j*l*iv.pi/N)**2/(N*q) for l,q in enumerate(p['P'],1))
        assert (spec-p['vv'][j-1]).a<=0<=(spec-p['vv'][j-1]).b
        print('variance spectral interval check',N,j,up(abs(spec-p['vv'][j-1])),flush=True)

if __name__=='__main__':
    for n in map(int,sys.argv[1:] or [80,320]):
        check(n);run(n);soft_run(n,'1e-5','.001')
