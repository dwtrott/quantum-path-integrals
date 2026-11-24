"""Interval checks of basis and classical inputs to projected_overlap.py (not BLAS)."""
from mpmath import iv
import mpmath, platform
import numpy as np
from projected_overlap import basis,geometry,target,ground
iv.dps=45
I=iv.mpf; Z=iv.mpc
J=Z(0,1)
def zz(z):return Z(float(np.real(z)),float(np.imag(z)))
def err(a,b):return abs(a-(zz(b) if isinstance(a,type(J)) else I(float(b)))).b
def conj(z):return Z(z.real,-z.imag)
def root(z):
    r=iv.sqrt((abs(z)+z.real)/2)
    assert r.a>0
    return Z(r,z.imag/(2*r))
def iground(y,v):return 1/(root(y)*root((1-J*v/y)/2))
def itarget(k,y,v):
    A=(1-J*v/y)/2
    z=(1/(A*abs(y)**2)-1)*(abs(y)/y)**2
    g=[Z(0) for _ in range(k)];g[0]=Z(1)
    for j in range(2,k,2):g[j]=g[j-2]*z*iv.sqrt(I(j-1)/j)
    a0=iv.sqrt(2/abs(v+J*y))
    tail=a0*abs(g[(k-1)//2*2])*abs(z)/iv.sqrt(1-abs(z)**2)
    return [a0*conj(x) for x in g],tail

def check_basis(k=48):
    x,a=basis(k);s=I(1)/16;norm=I(0)
    for i,xf in enumerate(x):
        t=I(float(xf));v=[iv.exp(-t*t/2)*iv.sqrt(s/iv.sqrt(iv.pi))]
        v.append(iv.sqrt(I(2))*t*v[0])
        for j in range(1,k-1):v.append(iv.sqrt(I(2)/(j+1))*t*v[j]-iv.sqrt(I(j)/(j+1))*v[j-1])
        for j in range(k):norm+=err(v[j],a[i,j])**2
    e=iv.sqrt(norm).b;print('basis Frobenius bound',e,flush=True);assert e< I('1e-11')

def check_geometry(n=1024,k=48,midpoint=True):
    y0=15/8 if midpoint else 1.;v0=(-3/16 if midpoint else 0.)+1j/y0
    rows,_,_,yf,vf=geometry(n,y0,v0)
    q,p=I(-2),I(2);y=Z(I(y0));v=Z(I(v0.real),1/I(y0));yi,vi=y,v
    d=I(2)/n;worst=I(0);rotmax=I(0)
    for j,(qf,yjf,wf,cf,mf) in enumerate(rows):
        m=d/2 if j in (0,n) else d
        w=3*iv.exp(-q*q/2);c=w*(1-q*q)
        assert abs(q).b<=I(2)+I('1e-20') and abs(y).b<(2 if midpoint else 3) and y.real.a>0
        unit=y/abs(y);uf=yjf/abs(yjf);r=Z(1);rf=1+0j
        for _ in range(k):
            rotmax=max(rotmax,err(r,rf));r*=unit;rf*=uf
        for h in (.05,.025,.0125):
            H=I(float(h));s=iv.sqrt(H)*abs(y);sf=np.sqrt(h)*abs(yjf)
            # h is defined as the binary64 input, not its rational decimal idealization.
            ef=err(q,qf);es=err(s,sf);ew=err(w,wf)
            ea=err(q*w*s,qf*wf*sf);eb=err(c*s*s/2,cf*sf*sf/2);eg=err(m/H,mf/h)
            assert s.b<1 and (m/H).b< I('.16')
            phase=(m/H).b*(2*(ef+24*es)+ew+24*ea+576*eb)+eg*2000
            worst=max(worst,phase.b)
        p+=m*q*w;v+=m*c*y
        if j<n:q+=d*p;y+=d*v
    print('geometry',midpoint,'parameter phase error',worst,'rotation max',rotmax,flush=True)
    assert worst<I('1e-9') and rotmax<I('1e-10')
    if midpoint:
        for K in (24,32,40,48):
            for yy,vv,yff,vff in ((yi,vi,y0,v0),(y,v,yf,vf)):
                g,t=itarget(K,yy,vv);gf,_,tf,_=target(K,yff,vff)
                e=iv.sqrt(sum(err(a,b)**2 for a,b in zip(g,gf))).b
                print('endpoint',K,'coefficient error',e,'tail error',err(t,tf),flush=True)
                assert e<I('1e-10') and err(t,tf)<I('1e-10')
    return y,v,yf,vf
if __name__=='__main__':
    print(platform.python_version(),np.__version__,mpmath.__version__,flush=True)
    check_basis()
    y,v,yf,vf=check_geometry();yo,vo,yof,vof=check_geometry(midpoint=False)
    yi=Z(I(15)/8);vi=Z(-I(3)/16,I(8)/15);a,b,c=iground(y,v),iground(yi,vi),iground(yo,vo)
    ratio=a*conj(b)/(c*abs(a*b));af,bf,cf=ground(yf,vf),ground(15/8,-3/16+1j/(15/8)),ground(yof,vof)
    rf=af*bf.conjugate()/(cf*abs(af*bf));print('ratio error',err(ratio,rf),'leading modulus error',err(abs(c),abs(cf)),flush=True)
    assert err(ratio,rf)<I('1e-10') and err(abs(c),abs(cf))<I('1e-10') and abs(c).a>I('.5')
