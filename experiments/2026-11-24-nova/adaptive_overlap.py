"""Standalone adaptive Hermite overlap test, based on Nova 051. IEEE conditional."""
import numpy as np
import time, platform

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
    lead=np.exp(1j*action/h)/(np.sqrt(y)*np.sqrt((1-1j*v/y)/2))
    return ref/lead

def run(n,h,ranks,delta=1e-8):
    _,_,_,yo,vo=geometry(n);y0=15/8;v0=-3/16+1j/y0
    rows,q,p,y,v=geometry(n,y0,v0);x,a=basis(48)
    # Contiguous matrices for each rank (avoid strided BLAS copies).
    mats={k:np.ascontiguousarray(a[:,:k]) for k in set(ranks)}
    b,tail=target(ranks[-1],y,v);f,tail0=target(ranks[0],y0,v0)
    ratio=ground(y,v)*ground(y0,v0).conjugate()/ground(yo,vo)
    ratio/=abs(ground(y,v)*ground(y0,v0));a0=abs(ground(yo,vo));fs=[];bs=[]
    def apply(c,row,sign,kout):
        assert np.linalg.norm(c)<2
        q,y,w,cc,m=row;s=np.sqrt(h)*abs(y);t=s*x
        r=3*np.exp(-(q+t)**2/2)-w+q*w*t+cc*t*t/2
        rot=np.empty(max(len(c),kout),complex);rot[0]=1;unit=y/abs(y)
        for j in range(1,len(rot)):rot[j]=rot[j-1]*unit
        wave=np.exp(-sign*1j*m/h*r)*(mats[len(c)]@(c*rot[:len(c)].conj()))
        u=rot[:kout]*(mats[kout].T@wave)
        resid=np.linalg.norm(wave-mats[kout]@(u*rot[:kout].conj()))
        return u,resid+delta
    for j,row in enumerate(rows):
        f,e=apply(f,row,1,ranks[j+1]);fs.append(e)
    val=ratio*np.vdot(b,f)
    for j in range(n,-1,-1):
        b,e=apply(b,rows[j],-1,ranks[j]);bs.append(e)
    bs=bs[::-1];prefix=tail0+delta;pair=0.
    for a,b in zip(fs,bs):
        pair+=prefix*b;prefix+=a
    pair+=prefix*(tail+delta)
    # No cumulative state error in residuals. Local projected orthogonality
    # defect <= delta, paired with backward vector of norm <2.
    radius=(pair+2*(n+3)*delta+10*delta)/(a0-delta)
    return val,radius*(1+1e-10)+10*delta,sum(fs),sum(bs)

def schedule(n,edge,mid,width):
    # symmetric high ranks at ends, low in middle; no pilots required
    return [edge if min(j,n+1-j)<=width else mid for j in range(n+2)]

if __name__=='__main__':
    print(platform.python_version(),np.__version__,flush=True)
    n=1024
    for h in (.05,.025,.0125):
        ref=fft(n,h);ref2=fft(n,h,32768);assert abs(ref-ref2)<1e-9
        for edge,mid,width in [(32,32,0),(40,40,0),(48,48,0),(48,24,64),(48,24,128),(48,24,256),(48,32,128),(48,32,256),(40,24,256),(40,32,128),(32,24,256)]:
            ranks=schedule(n,edge,mid,width);t=time.perf_counter();v,b,f,g=run(n,h,ranks);elapsed=time.perf_counter()-t
            assert abs(v-ref)<b
            print(h,edge,mid,width,'avg',sum(ranks)/len(ranks),'radius',b,'error',abs(v-ref),'seconds',elapsed,'f/g',f,g,flush=True)
