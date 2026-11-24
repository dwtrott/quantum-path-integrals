import mpmath as mp
from mpmath import iv
from fractions import Fraction as F
import math, sys, time
iv.dps=40
I=iv.mpf; J=iv.mpc(0,1)

def weights(n):
    # weights on [0,1], exact rational integration of Lagrange basis
    out=[]
    for j in range(n+1):
        p=[F(1)]
        for k in range(n+1):
            if k==j: continue
            z=[F(0)]*(len(p)+1)
            for i,c in enumerate(p):
                z[i]-=c*F(k,j-k); z[i+1]+=c*F(n,j-k)
            p=z
        out.append(sum(c/F(i+1) for i,c in enumerate(p)))
    assert all(sum(w*F(j,n)**k for j,w in enumerate(out))==F(1,k+1) for k in range(n+2))
    return out

def rat(x): return I(x.numerator)/x.denominator

def sinc(x):
    # valid for |x|<=1: omitted terms bounded by 2/(2p+3)!
    assert abs(x).b<=1
    p=15
    return sum((-1)**r*x**(2*r)/math.factorial(2*r+1) for r in range(p+1))+I([-1,1])*2/math.factorial(2*p+3)

def model(N,h,eta):
    alpha=iv.pi/N; delta=eta/N; theta=alpha+delta; mid=alpha+delta/2
    d=iv.sqrt(2)*iv.sin(theta/2); T=N*d
    a=-2*N/d*iv.sin(mid)*iv.sin(delta/2)
    assert abs(eta).b<=I('0.5')
    ds=sum((-1)**r*N*eta**(2*r-1)*(1-I(1)/(2*N)**(2*r))/math.factorial(2*r+1) for r in range(1,16))
    ds+=I([-1,1])*2*N*abs(eta)**31/math.factorial(33)
    K=1/(N*iv.sin(theta)*sinc(eta))
    D=(iv.cos(alpha+3*delta/4)*sinc(delta/4)*sinc(eta)/2+iv.sin(mid)*ds)/(N*iv.sin(mid)*sinc(delta/2)*iv.sin(theta)*sinc(eta))
    v=[iv.sin(j*alpha) for j in range(1,N)]
    P=[(h*d*v[j-1],-h*d*K*j*iv.cos(j*mid)*sinc(j*delta/2)) for j in range(1,N)]
    Q=[(D*v[j-1]-K*(N-j)*iv.cos((N-j)*mid)*sinc((N-j)*delta/2),iv.sin((N-j)*theta)) for j in range(1,N)]
    diag=[sum(p*q for p,q in zip(Pj,Qj)) for Pj,Qj in zip(P,Q)]
    # independent spectral check of diagonal at first, middle, last entries
    for j in set([1,N//2,N-1]):
        direct=sum(2*h/N*iv.sin(j*k*alpha)**2/((2/d)*(iv.cos(theta)-iv.cos(k*alpha))) for k in range(2,N))
        assert (diag[j-1]-direct).a<=0<=(diag[j-1]-direct).b
    tau=d/h; scale=iv.sqrt(iv.sqrt(h)); rot=iv.exp(-J*iv.pi/8); Qmax=I(1)/2
    def f(x):
        q=scale*rot*x; r=[q*y for y in v]; z=[-J*y*y/2 for y in r]; E=[iv.exp(y) for y in z]
        w=J*tau*sum(e-1-y for e,y in zip(E,z)); wp=[tau*y*(e-1) for y,e in zip(r,E)]; wpp=[tau*(e-1-J*y*y*e) for y,e in zip(r,E)]
        accum=[iv.mpc(0),iv.mpc(0)]; quad=iv.mpc(0)
        for pj,qj,xx in zip(P,Q,wp):
            quad+=xx*xx*sum(y*z for y,z in zip(pj,qj))+2*xx*sum(y*z for y,z in zip(accum,qj))
            accum=[y+xx*z for y,z in zip(accum,pj)]
        return rot*iv.exp(-a*q*q/(2*h)+w)*(1+(sum(x*y for x,y in zip(wpp,diag))+quad)/2)
    v4=sum(y**4 for y in v); M=max(x.b for x in diag)/h
    Gnorm=1/((2/d)*(iv.cos(theta)-iv.cos(2*alpha)))
    def disk(c,R):
        b=c+R; qmax=scale*b; assert qmax.b<=1
        lo=I(0) if (c-R).a<0 else c-R
        lower=lo**4-6*b*b*R*R
        expo=abs(a)*b*b/(2*iv.sqrt(h))-d*v4*lower/8+T*iv.sqrt(h)*b**6*iv.exp(qmax*qmax/2)/48
        return iv.exp(expo)*(1+(3*T*M*qmax*qmax+T*d*Gnorm*qmax**6/h)/2)
    return Qmax/scale,f,disk,T,a

def run(N,hs,es,panels,n=12):
    h=I(hs); eta=I(es); t=time.perf_counter(); end,f,disk,T,a=model(N,h,eta)
    W=weights(n); ww=[rat(x) for x in W]; l=end/(2*panels); R=8*l
    val=iv.mpc(0); err=I(0)
    for j in range(panels):
        c=(2*j+1)*l
        val+=2*l*sum(w*f(c+l*(I(2*k)/n-1)) for k,w in enumerate(ww))
        err+=2*l*(1+sum(rat(abs(w)) for w in W))*disk(c,R)*I(1)/8**(n+2)*I(8)/7
    val*=2; err*=2
    box=iv.mpc(val.real+I([-1,1])*err.b,val.imag+I([-1,1])*err.b)
    print('N,h,eta,panels',N,hs,es,panels,'T',T,'a',a,'value',val,'analytic radius',err,'enclosure',box,'sec',time.perf_counter()-t,flush=True)
    return box
if __name__=='__main__':
    print(sys.version,mp.__version__)
    for N,hs,es,panels in [(4,'1e-5','0',72),(4,'1e-5','0',96),(80,'1e-5','0.001',72),(80,'1e-5','-0.001',72)]: run(N,hs,es,panels)
