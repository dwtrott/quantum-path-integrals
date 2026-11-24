import math, time, sys
from fractions import Fraction as F
import mpmath as mp
from mpmath import iv
iv.dps=40
I=iv.mpf; J=iv.mpc(0,1)

def weights(n):
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
    assert abs(x).b<=1
    return sum((-1)**r*x**(2*r)/math.factorial(2*r+1) for r in range(16))+I([-1,1])*2/math.factorial(33)

def model(N,h,eta,p):
    alpha=iv.pi/N; delta=eta/N; theta=alpha+delta; mid=alpha+delta/2
    d=iv.sqrt(2)*iv.sin(theta/2); T=N*d; tau=d/h
    a=-2*N/d*iv.sin(mid)*iv.sin(delta/2)
    assert abs(eta).b<=I('0.5')
    ds=sum((-1)**r*N*eta**(2*r-1)*(1-I(1)/(2*N)**(2*r))/math.factorial(2*r+1) for r in range(1,16))
    ds+=I([-1,1])*2*N*abs(eta)**31/math.factorial(33)
    K=1/(N*iv.sin(theta)*sinc(eta))
    D=(iv.cos(alpha+3*delta/4)*sinc(delta/4)*sinc(eta)/2+iv.sin(mid)*ds)/(N*iv.sin(mid)*sinc(delta/2)*iv.sin(theta)*sinc(eta))
    v=[iv.sin(j*alpha) for j in range(1,N)]
    P=[(h*d*v[j-1],-h*d*K*j*iv.cos(j*mid)*sinc(j*delta/2)) for j in range(1,N)]
    Q=[(D*v[j-1]-K*(N-j)*iv.cos((N-j)*mid)*sinc((N-j)*delta/2),iv.sin((N-j)*theta)) for j in range(1,N)]
    diag=[sum(x*y for x,y in zip(u,w)) for u,w in zip(P,Q)]
    Gnorm=1/((2/d)*(iv.cos(theta)-iv.cos(2*alpha)))
    for j in set([1,N//2,N-1]):
        direct=sum(2*h/N*iv.sin(j*k*alpha)**2/((2/d)*(iv.cos(theta)-iv.cos(k*alpha))) for k in range(2,N))
        assert (diag[j-1]-direct).a<=0<=(diag[j-1]-direct).b
    c=[(-J/2)**k/math.factorial(k) for k in range(p+1)]
    # powers shared across the polynomial coefficients
    ev=[[x**(2*k) for k in range(1,p+1)] for x in v]
    od=[[x*y for y in row] for x,row in zip(v,ev)]
    wc=[iv.mpc(0)]*(p+1); bc=[iv.mpc(0)]*(2*p+2)
    for k in range(1,p+1):
        if k>=2: wc[k]=J*tau*c[k]*sum(row[k-1] for row in ev)
        bc[k]=tau*(2*k+1)*c[k]*sum(row[k-1]*g for row,g in zip(ev,diag))/2
    # Accumulate all real bilinear contractions of odd powers with hG, using its rank-two generators.
    # diagonal plus both triangles; combine polynomial coefficients by k+l.
    acc=[[I(0),I(0)] for _ in range(p)]
    M=[[I(0) for _ in range(p)] for _ in range(p)]
    for row,pj,qj,g in zip(od,P,Q,diag):
        b=[sum(x*y for x,y in zip(z,qj)) for z in acc]
        for k in range(p):
            for l in range(k,p):
                M[k][l]+=g*row[k]*row[l]+b[k]*row[l]+b[l]*row[k]
        for k in range(p):
            acc[k]=[x+row[k]*y for x,y in zip(acc[k],pj)]
    for k in range(1,p+1):
        for l in range(k,p+1):
            bc[k+l+1]+=tau*tau*c[k]*c[l]*M[k-1][l-1]*(1 if k==l else 2)/2
    bc[0]=iv.mpc(1)
    scale=iv.sqrt(iv.sqrt(h)); rot=iv.exp(-J*iv.pi/8)
    def horner(C,z):
        y=C[-1]
        for x in reversed(C[:-1]): y=y*z+x
        return y
    def f(x):
        q=scale*rot*x; z=q*q
        return rot*iv.exp(-a*z/(2*h)+horner(wc,z))*horner(bc,z)
    # Pointwise bound on exact f minus polynomial-exponential f on the real scaled axis.
    def point_error(x):
        q=scale*x; b=q*q/2
        R=iv.exp(b)*b**(p+1)/math.factorial(p+1)
        ew=T/h*R
        # w' error <= tau |q| R per coordinate; |w'| <= tau |q| b exp(b).
        e1=tau*q*R; w1=tau*q*b*iv.exp(b)
        # w''=tau[(exp(z)-1)+2z exp(z)], truncated with (2k+1)c_k.
        e2=tau*(R+2*b*iv.exp(b)*b**p/math.factorial(p))
        w2=tau*3*b*iv.exp(b)
        tr=sum(diag)
        # |v^T C v| <= ||C|| ||v||_2^2 also for complex bilinear v.
        eb=(tr*e2+h*Gnorm*(N-1)*e1*(2*w1+e1))/2
        B=1+(tr*w2+h*Gnorm*(N-1)*w1*w1)/2
        # From Nova 029, on the ray |q|<=1, Re w <= -3T q^4/(128h).
        mag=iv.exp(-a*q*q/(2*iv.sqrt(2)*h)-3*T*q**4/(128*h))
        return mag*iv.exp(ew)*(eb+B*ew)
    v4=sum(y**4 for y in v); Mdiag=max(x.b for x in diag)/h
    def disk(c,R):
        b=c+R; qmax=scale*b; assert qmax.b<=1
        lo=I(0) if (c-R).a<0 else c-R
        lower=lo**4-6*b*b*R*R
        expo=abs(a)*b*b/(2*iv.sqrt(h))-d*v4*lower/8+T*iv.sqrt(h)*b**6*iv.exp(qmax*qmax/2)/48
        return iv.exp(expo)*(1+(3*T*Mdiag*qmax*qmax+T*d*Gnorm*qmax**6/h)/2)
    return f,point_error,disk,T,a,Gnorm,[x/h for x in diag],v

def run(N,hs,es,panels,p=6,n=12,Qstr='0.3'):
    h=I(hs); eta=I(es); t=time.perf_counter(); f,pe,disk,T,a,Gnorm,G,v=model(N,h,eta,p)
    setup=time.perf_counter()-t
    scale=iv.sqrt(iv.sqrt(h)); Q=I(Qstr); end=Q/scale
    W=weights(n); ww=[rat(x) for x in W]; aw=[rat(abs(x)) for x in W]; l=end/(2*panels); R=8*l
    val=iv.mpc(0); err=I(0); ep=I(0)
    for j in range(panels):
        c=(2*j+1)*l; nodes=[c+l*(I(2*k)/n-1) for k in range(n+1)]
        val+=2*l*sum(w*f(x) for w,x in zip(ww,nodes))
        ep+=2*l*sum(w*pe(abs(x)) for w,x in zip(aw,nodes))
        err+=2*l*(1+sum(aw))*disk(c,R)*I(1)/8**(n+2)*I(8)/7
    val*=2; err*=2; ep*=2
    d=T/N; c=3*T/128; assert abs(a).b <= (4*c*Q**2).a
    tail=2*(I('0.5')-Q)/scale*iv.exp(abs(a)*Q**2/(2*h)-c*Q**4/h)*(1+(3*T*(3*T/4)+T*d*Gnorm/h)/2)
    print('N,h,eta,panels,p',N,hs,es,panels,p,'value',val,'quadrature',err,'polynomial',ep,'tail',tail,'setup',setup,'seconds',time.perf_counter()-t,flush=True)
    return val,err,ep,tail,(T,a,Gnorm,G,v)
if __name__=='__main__':
    print(sys.version,mp.__version__)
    for N,e,panels in [(4,'0',24),(80,'-0.001',24),(80,'0.001',18),(80,'0.001',24),(320,'0.001',24)]: run(N,'1e-5',e,panels)
