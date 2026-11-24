"""Projected real-target kicks. Run OPENBLAS_NUM_THREADS=1 python projected_overlap.py."""
import math, time, sys, platform
import numpy as np


def geometry(n,y0=1+0j,v0=1j):
    d=2./n; q=-2.; p=2.; y=y0; v=v0; rows=[]
    for j in range(n+1):
        m=d/2 if j in (0,n) else d
        w=3*np.exp(-q*q/2); c=w*(1-q*q)
        rows.append((q,y,w,c,m))
        p+=m*q*w; v+=m*c*y
        if j<n: q+=d*p; y+=d*v
    return rows,q,p,y,v


def basis(k,step=1/16,B=24):
    x=np.arange(-round(B/step),round(B/step)+1)*step
    a=np.empty((len(x),k)); a[:,0]=np.exp(-x*x/2)*np.sqrt(step/np.sqrt(np.pi))
    a[:,1]=np.sqrt(2)*x*a[:,0]
    for j in range(1,k-1):a[:,j+1]=np.sqrt(2/(j+1))*x*a[:,j]-np.sqrt(j/(j+1))*a[:,j-1]
    return x,a


def target(k,y,v):
    A=(1-1j*v/y)/2
    z=(1/(A*abs(y)**2)-1)*(abs(y)/y)**2
    g=np.zeros(k,complex);g[0]=1
    for j in range(2,k,2):g[j]=g[j-2]*z*np.sqrt((j-1)/j)
    a0=np.sqrt(2/abs(v+1j*y))
    # |g| has been scaled to unit full norm (and common phase removed).
    tail=a0*abs(g[(k-1)//2*2])*abs(z)/np.sqrt(1-abs(z)**2)
    return (a0*g).conj(),a0,tail,abs(z)


def ground(y,v):return 1/(np.sqrt(y)*np.sqrt((1-1j*v/y)/2))


def propagate(n,h,k,step=1/16,delta=1e-8,midpoint=True):
    _,_,_,yold,vold=geometry(n)
    # Fixed rounded midpoint, exact real parameters.
    y0=15/8 if midpoint else 1.;v0=(-3/16 if midpoint else 0.)+1j/y0
    rows,q,p,y,v=geometry(n,y0,v0);x,a=basis(k,step)
    b,_,tail,z=target(k,y,v)
    f,_,tail0,_=target(k,y0,v0)
    tail+=delta;tail0+=delta
    ratio=ground(y,v)*ground(y0,v0).conjugate()/ground(yold,vold)
    ratio/=abs(ground(y,v)*ground(y0,v0))
    a0=abs(ground(yold,vold));fs=[];bs=[]
    def local(row,sign):
        q,y,w,c,m=row;s=np.sqrt(h)*abs(y);t=s*x
        # R(q+t) = V(q+t)-[V(q)+V'(q)t+V''(q)t^2/2]
        r=3*np.exp(-(q+t)**2/2)-w+q*w*t+c*t*t/2
        rot=np.empty(k,complex);rot[0]=1;unit=y/abs(y)
        for j in range(1,k):rot[j]=rot[j-1]*unit
        return np.exp(-sign*1j*m/h*r),rot
    def apply(c,row,sign):
        assert np.linalg.norm(c)<2
        phase,rot=local(row,sign)
        wave=phase*(a@(c*rot.conj()))
        u=rot*(a.T@wave)
        # Conditional direct residual allowance: avoid subtracting two nearly equal norms.
        resid=np.linalg.norm(wave-a@(u*rot.conj()))
        return u,resid+delta
    for row in rows:
        f,e=apply(f,row,1);fs.append(e)
    val=ratio*np.vdot(b,f)
    for row in rows[::-1]:
        b,e=apply(b,row,-1);bs.append(e)
    bs=bs[::-1];prefix=pair=0.
    for j,(a,b) in enumerate(zip(fs,bs)):
        # Difference between exact and computed projected state <= number of prior kicks times delta.
        pair+=prefix*(b+(n-j+1)*delta)
        prefix+=a+(j+1)*delta
    bound=(pair+tail*prefix+tail0*(sum(bs)+(n+1)*(n+2)/2*delta)+tail*tail0+(n+5)*delta)/(a0-delta)
    return val,bound*(1+1e-10)+10*delta,prefix,sum(bs),tail,z


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
    lead=np.exp(1j*action/h)/(np.sqrt(y)*np.sqrt((1-1j*v/y)/2))
    return ref/lead

if __name__=='__main__':
    print(platform.python_version(),np.__version__,flush=True)
    for h in (.05,.025,.0125):
        n=1024
        t=time.perf_counter();ref=fft(n,h);tf=time.perf_counter()-t
        refined=fft(n,h,32768)
        assert abs(ref-refined)<1e-9
        for k in (24,32,40,48):
            t=time.perf_counter();val,bound,f,b,tail,z=propagate(n,h,k);tp=time.perf_counter()-t
            print(n,h,k,'error',abs(val-ref),'conditional_bound',bound,'forward',f,'backward',b,'tail',tail,'seconds',tp,'fft_seconds',tf,'fft_change',abs(ref-refined),flush=True)
            assert abs(val-ref)<bound+1e-9
        t=time.perf_counter();val,bound,*_=propagate(n,h,32,midpoint=False)
        print('unsqueezed',h,abs(val-ref),bound,'seconds',time.perf_counter()-t,flush=True)
