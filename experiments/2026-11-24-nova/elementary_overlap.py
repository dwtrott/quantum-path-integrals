"""Nova 053: explicit elementary kernels and no-BLAS residual propagation.
Run OPENBLAS_NUM_THREADS=1 python elementary_overlap.py (assertions on).
The FFT comparison is floating point, not part of the certificate.
"""
import math, platform, time
import numpy as np
from mpmath import iv
import mpmath

# All polynomial coefficients rounded from exact integer factorials.
EXP = [(-1.)**j / math.factorial(j) for j in range(17)]
COS = [(-1.)**j / math.factorial(2*j) for j in range(11)]
SIN = [(-1.)**j / math.factorial(2*j+1) for j in range(10)]

def horner(x, coefficients):
    v = np.full_like(x, coefficients[-1], dtype=np.float64)
    for c in coefficients[-2::-1]:
        v = v*x+c
    return v

def expneg(t):
    assert np.all((t >= 0) & (t <= 400))
    v = horner(np.minimum(t,40.)/64., EXP)
    for _ in range(6):
        v = v*v
    return v

def cis(t):
    assert np.all(np.abs(t) <= 150)
    q = np.rint(t/(np.pi/2))
    s = t-q*(np.pi/2)
    assert np.all(np.abs(s) < .786)
    c = horner(s*s,COS); a = s*horner(s*s,SIN)
    k = q.astype(np.int64)%4
    return np.where(k==0,c,np.where(k==1,-a,np.where(k==2,-c,a))) + 1j*np.where(k==0,a,np.where(k==1,c,np.where(k==2,-a,-c)))

def norm(x):
    return np.sqrt(np.sum(x.real*x.real+x.imag*x.imag))

def geometry(n,y0=1+0j,v0=1j):
    d=2./n; q=-2.; p=2.; y=y0; v=v0; rows=[]
    for j in range(n+1):
        m=d/2 if j in (0,n) else d
        w=3*np.exp(-q*q/2); c=w*(1-q*q)
        rows.append((q,y,w,c,m))
        p+=m*q*w; v+=m*c*y
        if j<n: q+=d*p; y+=d*v
    return rows,q,p,y,v

def basis(k):
    step=1/16; x=np.arange(-384,385)*step
    a=np.empty((len(x),k)); a[:,0]=np.exp(-x*x/2)*np.sqrt(step/np.sqrt(np.pi))
    a[:,1]=np.sqrt(2)*x*a[:,0]
    for j in range(1,k-1):a[:,j+1]=np.sqrt(2/(j+1))*x*a[:,j]-np.sqrt(j/(j+1))*a[:,j-1]
    return x,a

def target(k,y,v):
    A=(1-1j*v/y)/2; z=(1/(A*abs(y)**2)-1)*(abs(y)/y)**2
    g=np.zeros(k,complex);g[0]=1
    for j in range(2,k,2):g[j]=g[j-2]*z*np.sqrt((j-1)/j)
    a0=np.sqrt(2/abs(v+1j*y));tail=a0*abs(g[(k-1)//2*2])*abs(z)/np.sqrt(1-abs(z)**2)
    return (a0*g).conj(),tail

def ground(y,v):return 1/(np.sqrt(y)*np.sqrt((1-1j*v/y)/2))

def fft(n,h,M=16384):
    rows,q,p,y,v=geometry(n);d=2./n;L=40.;dx=2*L/M;x=np.arange(M)*dx-L
    k=2*np.pi*np.fft.fftfreq(M,dx);w=3*np.exp(-x*x/2)
    psi=np.exp(-(x+2)**2/(2*h)+2j*x/h)/(np.pi*h)**.25
    half=np.exp(-.5j*d*w/h);drift=np.exp(-.5j*h*d*k*k)
    psi*=half
    for j in range(n):
        psi=np.fft.ifft(drift*np.fft.fft(psi));psi*=half if j==n-1 else half*half
    ref=dx*np.dot(np.exp(-(x-q)**2/(2*h)-1j*p*x/h)/(np.pi*h)**.25,psi)
    action=-4-p*q+sum((b[0]-a[0])**2/(2*d) for a,b in zip(rows,rows[1:]))-sum(r[2]*r[4] for r in rows)
    lead=np.exp(1j*action/h)*ground(y,v)
    return ref/lead

iv.dps=45
I=iv.mpf; Z=iv.mpc; J=Z(0,1)
def zz(z):return Z(float(np.real(z)),float(np.imag(z)))
def err(a,b):return abs(a-(zz(b) if isinstance(a,type(J)) else I(float(b)))).b
def conj(z):return Z(z.real,-z.imag)
def root(z):
    r=iv.sqrt((abs(z)+z.real)/2); assert r.a>0
    return Z(r,z.imag/(2*r))
def iground(y,v):return 1/(root(y)*root((1-J*v/y)/2))
def itarget(k,y,v):
    A=(1-J*v/y)/2; z=(1/(A*abs(y)**2)-1)*(abs(y)/y)**2
    g=[Z(0) for _ in range(k)];g[0]=Z(1)
    for j in range(2,k,2):g[j]=g[j-2]*z*iv.sqrt(I(j-1)/j)
    a0=iv.sqrt(2/abs(v+J*y))
    tail=a0*abs(g[(k-1)//2*2])*abs(z)/iv.sqrt(1-abs(z)**2)
    return [a0*conj(x) for x in g],tail

def validate(n,ks):
    x,a=basis(max(ks));s=I(1)/16;total=I(0)
    for i,xf in enumerate(x):
        t=I(float(xf));v=[iv.exp(-t*t/2)*iv.sqrt(s/iv.sqrt(iv.pi))]
        v.append(iv.sqrt(I(2))*t*v[0])
        for j in range(1,max(ks)-1):v.append(iv.sqrt(I(2)/(j+1))*t*v[j]-iv.sqrt(I(j)/(j+1))*v[j-1])
        for j in range(max(ks)):total+=err(v[j],a[i,j])**2
    e=iv.sqrt(total).b;print('basis error',e,flush=True);assert e<I('1e-11')
    for k in ks:assert np.array_equal(basis(k)[1],a[:,:k])
    gs=[];gfs=[]
    for mid in (True,False):
        y0=15/8 if mid else 1.;v0=(-3/16 if mid else 0.)+1j/y0
        rows,_,_,yf,vf=geometry(n,y0,v0)
        q,p=I(-2),I(2);y=Z(I(y0));v=Z(I(v0.real),1/I(y0));yi,vi=y,v
        d=I(2)/n;worst=I(0);rotmax=I(0)
        for j,(qf,yjf,wf,cf,mf) in enumerate(rows):
            m=d/2 if j in (0,n) else d
            w=3*iv.exp(-q*q/2);c=w*(1-q*q)
            assert abs(q).b<=I(2)+I('1e-20') and abs(y).b<(2 if mid else 3) and y.real.a>0
            unit=y/abs(y);uf=yjf/abs(yjf);r=Z(1);rf=1+0j
            for _ in range(max(ks)):
                rotmax=max(rotmax,err(r,rf));r*=unit;rf*=uf
            for h in (.05,.025,.0125):
                H=I(float(h));ss=iv.sqrt(H)*abs(y);sf=np.sqrt(h)*abs(yjf)
                ef=err(q,qf);es=err(ss,sf);ew=err(w,wf)
                ea=err(q*w*ss,qf*wf*sf);eb=err(c*ss*ss/2,cf*sf*sf/2);eg=err(m/H,mf/h)
                assert ss.b<1 and (m/H).b<I('.16')
                phase=(m/H).b*(2*(ef+24*es)+ew+24*ea+576*eb)+eg*2000
                worst=max(worst,phase.b)
            p+=m*q*w;v+=m*c*y
            if j<n:q+=d*p;y+=d*v
        print('geometry',mid,'phase',worst,'rotation',rotmax,flush=True)
        assert worst<I('1e-10') and rotmax<I('1e-11')
        if mid:
            for k in ks:
                for yy,vv,yff,vff in ((yi,vi,y0,v0),(y,v,yf,vf)):
                    g,t=itarget(k,yy,vv);gf,tf=target(k,yff,vff)
                    e=iv.sqrt(sum(err(a,b)**2 for a,b in zip(g,gf))).b
                    assert e<I('1e-10') and err(t,tf)<I('1e-10')
            gs.extend((iground(y,v),iground(yi,vi)));gfs.extend((ground(yf,vf),ground(y0,v0)))
        else:
            gs.append(iground(y,v));gfs.append(ground(yf,vf))
    a,b,c=gs;af,bf,cf=gfs
    ratio=a*conj(b)/(c*abs(a*b));rf=af*bf.conjugate()/(cf*abs(af*bf))
    print('ratio error',err(ratio,rf),'leading modulus error',err(abs(c),abs(cf)),flush=True)
    assert err(ratio,rf)<I('1e-10') and err(abs(c),abs(cf))<I('1e-10') and abs(c).a>I('.5')
    # Interval sanity checks of polynomial constants and analytic quadrature bound.
    for coeff,degree in ((EXP,lambda j:j),(COS,lambda j:2*j),(SIN,lambda j:2*j+1)):
        for j,c in enumerate(coeff):assert err(I((-1)**j)/math.factorial(degree(j)),c)<I('1e-16')
    assert err(iv.pi,np.pi)<I('2.23e-16')
    assert (2*iv.exp(24)/(iv.exp(32*iv.pi)-1)+iv.exp(-230)).b<I('1e-31')

# No transcendental or BLAS calls inside explicit propagation. Reductions
# sum at most 769 terms; the error bound permits any binary addition tree.
def matvec(a,c):
    return np.sum(a*c[None,:],axis=1)
def tmatvec(a,c):
    return np.sum(a*c[:,None],axis=0)

def run(n,h,k,mode='explicit',delta=1e-8):
    _,_,_,yo,vo=geometry(n);y0=15/8;v0=-3/16+1j/y0
    rows,q,p,y,v=geometry(n,y0,v0);x,a=basis(k)
    b,tail=target(k,y,v);f,tail0=target(k,y0,v0)
    ratio=ground(y,v)*ground(y0,v0).conjugate()/ground(yo,vo)
    ratio/=abs(ground(y,v)*ground(y0,v0));a0=abs(ground(yo,vo));fs=[];bs=[]
    # Precompute identical phase and rotation for forward and backward passes.
    phases=[];rots=[]
    for q,y,w,cc,m in rows:
        s=np.sqrt(h)*abs(y);t=s*x;arg=(q+t)**2/2
        r=3*(np.exp(-arg) if mode=='libm' else expneg(arg))-w+q*w*t+cc*t*t/2
        angle=-(m/h)*r
        phases.append(np.exp(1j*angle) if mode=='libm' else cis(angle))
        rot=np.empty(k,complex);rot[0]=1;unit=y/abs(y)
        for j in range(1,k):rot[j]=rot[j-1]*unit
        rots.append(rot)
    def apply(c,j,sign):
        assert norm(c)<2
        rot=rots[j];phase=phases[j] if sign==1 else phases[j].conj()
        if mode=='explicit':
            wave=phase*matvec(a,c*rot.conj());u=rot*tmatvec(a,wave)
            resid=norm(wave-matvec(a,u*rot.conj()))
        else:
            wave=phase*(a@(c*rot.conj()));u=rot*(a.T@wave)
            resid=norm(wave-a@(u*rot.conj()))
        return u,resid+delta
    for j in range(n+1):
        f,e=apply(f,j,1);fs.append(e)
    assert norm(f)<2 and norm(b)<2
    val=ratio*np.sum(b.conj()*f)
    for j in range(n,-1,-1):
        b,e=apply(b,j,-1);bs.append(e)
    assert norm(b)<2
    prefix=tail0+delta;pair=0.
    for r,s in zip(fs,bs[::-1]):pair+=prefix*s;prefix+=r
    pair+=prefix*(tail+delta)
    radius=(pair+2*(n+3)*delta+10*delta)/(a0-delta)
    return val,radius*(1+1e-10)+10*delta

if __name__=='__main__':
    print(platform.python_version(),np.__version__,mpmath.__version__,flush=True)
    n=1024;ks=(32,40,48)
    t=time.perf_counter();validate(n,ks);print('validation seconds',time.perf_counter()-t,flush=True)
    # Independent interval checks on a modest adversarial grid; not the uniform proof.
    t=np.r_[np.linspace(0,400,1001),40-np.spacing(40),40+np.spacing(40)]
    e=expneg(t);ee=max(err(iv.exp(-I(float(x))),y) for x,y in zip(t,e))
    t=np.r_[np.linspace(-150,150,1001),np.arange(-94,95)*np.pi/2+np.pi/4]
    c=cis(t);ce=max(err(Z(iv.cos(I(float(x))),iv.sin(I(float(x)))),y) for x,y in zip(t,c))
    print('interval sample exp/cis',ee,ce,flush=True);assert ee<I('1e-12') and ce<I('1e-12')
    for h in (.05,.025,.0125):
        t=time.perf_counter();ref=fft(n,h);ft=time.perf_counter()-t
        ref2=fft(n,h,32768);assert abs(ref-ref2)<1e-9
        for k in ks:
            vals=[]
            for mode in ('libm','taylor','explicit'):
                t=time.perf_counter();v,b=run(n,h,k,mode);elapsed=time.perf_counter()-t
                assert abs(v-ref)<b;vals.append(v)
                print(h,k,mode,'radius',b,'error',abs(v-ref),'seconds',elapsed,'fft seconds',ft,flush=True)
            assert max(abs(v-vals[0]) for v in vals)<1e-9
