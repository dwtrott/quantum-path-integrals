"""Linear Green contractions for Sparks 048, real classical target only.
Run: OPENBLAS_NUM_THREADS=1 python green_a1.py
"""
import sys, time, platform
import numpy as np
import mpmath as mp
from mpmath import iv


def linear(n, interval=False):
    # The same recurrence works in complex arithmetic or complex intervals.
    C, exp = (iv.mpc, iv.exp) if interval else (complex, np.exp)
    d = iv.mpf(2)/n if interval else 2./n
    m = [d/2] + [d]*(n-1) + [d/2]
    x = [-2]; p = 2 + 3*m[0]*x[0]*exp(-x[0]**2/2)
    for j in range(n):
        x.append(x[-1]+d*p)
        p += 3*m[j+1]*x[-1]*exp(-x[-1]**2/2)
    w = [exp(-v*v/2) for v in x]
    c = [3*(1-v*v)*a for v,a in zip(x,w)]
    def shoot(c):
        y, p = C(1), C(0,1)
        out = [y]
        for c0,c1 in zip(c,c[1:]):
            p += d*c0*y/2
            y += d*p
            p += d*c1*y/2
            out.append(y)
        return out, p+C(0,1)*y
    l, W = shoot(c)
    u, _ = shoot(c[::-1]); u.reverse()
    if interval:
        assert float((W.real**2+W.imag**2).a)>0
    gd = [a*b/W for a,b in zip(l,u)]
    w3 = [a*(3*v-v**3)*b for a,v,b in zip(m,x,w)]
    w4 = [a*(v**4-6*v*v+3)*b for a,v,b in zip(m,x,w)]
    b = [v*z for v,z in zip(w3,gd)]
    # Strict prefixes: each pair i<j contributes twice, not conjugated.
    t1=t3=s1=s3=C(0)
    for a,z,bj,cj in zip(l,u,b,w3):
        t1 += 2*bj*z*s1/W + bj*bj*a*z/W
        t3 += 2*cj*z**3*s3/W**3 + cj*cj*(a*z/W)**3
        s1 += bj*a
        s3 += cj*a**3
    a1 = C(0,3)/8*sum(v*z*z for v,z in zip(w4,gd)) + C(0,9)/72*(9*t1+6*t3)
    # Two-sided unitary-product bound; no absolute contour determinant loss.
    # Use M4=9 and M3^2<23 for V=3 exp(-x^2/2).
    sq = (lambda z: z.real**2+z.imag**2) if interval else (lambda z: abs(z)**2)
    root = iv.sqrt if interval else np.sqrt
    S4=S6=prefix=0
    for mj,lj,uj in zip(m,l,u):
        L,U=sq(lj),sq(uj)
        S4 += mj*L*U
        f,bk=mj*L*root(L),mj*U*root(U)
        S6 += bk*(f+2*prefix)
        prefix += f
    bound = (9*S4/32+15*23*S6/576)*root(root(sq(W)))/root(2)
    return a1, W, x, m, bound, p



def dense(n):
    # Independent construction of the kinetic matrix and full inverse.
    d=2./n; m=np.full(n+1,d); m[[0,-1]]/=2
    x=np.empty(n+1); x[0]=-2.; p=2.
    for j in range(n):
        p += .5*d*3*x[j]*np.exp(-x[j]**2/2)
        x[j+1]=x[j]+d*p
        p += .5*d*3*x[j+1]*np.exp(-x[j+1]**2/2)
    w=np.exp(-x*x/2)
    H=(np.diag(np.r_[1,2*np.ones(n-1),1])-np.diag(np.ones(n),1)-np.diag(np.ones(n),-1))/d
    H=H.astype(complex)-np.diag(3*m*(x*x-1)*w)
    H[0,0]+=1j; H[-1,-1]+=1j
    G=np.linalg.inv(H); gd=np.diag(G)
    b=m*(3*x-x**3)*w; c=m*(x**4-6*x*x+3)*w
    a=3j/8*(c@gd**2)+9j/72*(9*(b*gd)@G@(b*gd)+6*b@(G**3)@b)
    sign, logdet=np.linalg.slogdet(H)
    W=sign*np.exp(logdet+n*np.log(d))
    return a,W


def leading(n,h,x,m,p):
    # Exact quadratic propagation, no determinant sign fitting.
    d=2./n; y=1+0j; v=1j; phase=0.
    for j in range(n):
        cj=3*(1-x[j]**2)*np.exp(-x[j]**2/2)
        v += m[j]*cj*y
        yn=y+d*v
        phase += np.angle(yn/y)
        y=yn
    v += m[-1]*3*(1-x[-1]**2)*np.exp(-x[-1]**2/2)*y
    W=v+1j*y
    phase += np.angle(W/y)
    phi=2*x[0]-p*x[-1]+sum((b-a)**2/(2*d) for a,b in zip(x,x[1:]))-3*sum(mj*np.exp(-v*v/2) for v,mj in zip(x,m))
    return np.sqrt(2/abs(W))*np.exp(1j*(phi/h+np.pi/4-phase/2))


def fft_reference(n,h,q,p,M):
    # Normalized coherent-state overlap; same time slicing, independent spatial grid.
    d=2./n; L=40.; dx=2*L/M; x=np.arange(M)*dx-L
    k=2*np.pi*np.fft.fftfreq(M,dx)
    psi=np.exp(-(x+2)**2/(2*h)+2j*x/h)/(np.pi*h)**.25
    half=np.exp(-.5j*d*3*np.exp(-x*x/2)/h)
    drift=np.exp(-.5j*h*d*k*k)
    psi=half*psi
    for j in range(n):
        psi=np.fft.ifft(drift*np.fft.fft(psi))
        psi*=half if j==n-1 else half*half
    return dx*np.dot(np.exp(-(x-q)**2/(2*h)-1j*p*x/h)/(np.pi*h)**.25,psi)


def median_time(f,n,k=3):
    times=[]
    for _ in range(k):
        t=time.perf_counter(); f(n); times.append(time.perf_counter()-t)
    return float(np.median(times))


if __name__=='__main__':
    print(platform.python_version(),np.__version__,mp.__version__)
    iv.dps=50
    for n in (1,2,16,64,256,1024):
        a,W,*_=linear(n); ref,V=dense(n)
        assert abs(a-ref)<2e-8 and abs(W-V)<2e-7,(n,a,ref,W,V)
        t=median_time(linear,n); q=median_time(dense,n)
        print('dense_check',n,'error',abs(a-ref),'W_error',abs(W-V),'linear_s',t,'dense_s',q,flush=True)
    for n in (16,64,256,1024,4096):
        t=time.perf_counter(); a,W,_,_,bound,_=linear(n,True)
        # Print interval endpoints (the full values, not a float center/radius).
        print('interval',n,'a1',a,'W',W,'relative_C',bound,'seconds',time.perf_counter()-t,flush=True)
        assert float((a.real.delta+a.imag.delta).b)<1e-25
    for n in (4096,16384,65536):
        a,W,*_=linear(n)
        print('large_float',n,a,'W',W,'linear_s',median_time(linear,n),flush=True)

    for n in (16,64,256,1024):
        a,W,x,m,bound,p=linear(n)
        for h in (.05,.025,.0125):
            t=time.perf_counter(); L=leading(n,h,x,m,p); dt=time.perf_counter()-t
            t=time.perf_counter(); ref=fft_reference(n,h,x[-1],p,16384); ft=time.perf_counter()-t
            refined=fft_reference(n,h,x[-1],p,32768)
            assert abs(ref-refined)<2e-8 and abs(ref/L-1)<h*bound
            print('fft',n,h,'C',bound,'error',abs(ref/L-1),'corrected_error',abs(ref/L-1-h*a),'grid_change',abs(ref-refined),'leading_s',dt,'fft_s',ft,flush=True)
