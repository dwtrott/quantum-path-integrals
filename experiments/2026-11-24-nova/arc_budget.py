import math, sys, mpmath as mp
from mpmath import iv
iv.dps=40
I=iv.mpf

def upper(x): return float(x.b)
def polyadd(a,b):
    c=a.copy()
    for k,v in b.items(): c[k]=c.get(k,I(0))+v
    return c

def polymul(a,b):
    c={}
    for (i,j),v in a.items():
        for (k,l),w in b.items(): c[i+k,j+l]=c.get((i+k,j+l),I(0))+v*w
    return c

def scale(a,c): return {k:v*c for k,v in a.items()}
def gm(k):
    # odd moments bounded by Cauchy-Schwarz; even moments exact.
    if k%2==0: return I(math.prod(range(1,k,2)))
    return iv.sqrt(gm(k-1)*gm(k+1))

def run(N,hs,es,variant="sharp"):
    h=I(hs); eta=I(es); alpha=iv.pi/N; theta=(iv.pi+eta)/N; d=iv.sqrt(2)*iv.sin(theta/2); T=N*d
    a=-2*N/d*iv.sin(alpha+eta/(2*N))*iv.sin(eta/(2*N)); k=abs(a)/iv.sqrt(h)
    v=[iv.sin(j*alpha) for j in range(1,N)]
    lam=[(2/d)*(iv.cos(theta)-iv.cos(j*alpha)) for j in range(2,N)]
    assert min(x.a for x in lam)>0 and abs(eta).b<=I('0.5') and N>=4
    # Use a diagonal-only independent spectral computation, O(N^2); not production.
    G=[sum(2*iv.sin(j*l*alpha)**2/(N*L) for l,L in zip(range(2,N),lam)) for j in range(1,N)]
    M=max(x.b for x in G); g=1/lam[0]; delta=I(1)/10 if variant=="baseline" else I(35)/1000; K=20 if variant=="baseline" else 626; c4=I(1)/(32 if variant=="baseline" else 16); s=iv.sqrt(iv.sqrt(h))
    A=T/512; B=T/128
    assert abs(a).b<=min(A.a,B.a)
    # Tail bound by tangent lower bounds on q^6 and q^2.
    ct=A/(2*h); cg=B/(2*h); rt=iv.sqrt(2)
    tail=2*iv.exp(-ct)/(6*ct)+2*iv.exp(-2*cg)/(2*cg*rt)
    # W from 027, with Gamma(1/6)/3 < 2 (verified by elementary integral below).
    W=iv.exp(k**I('1.5')*h**I('-.25')/(3*iv.sqrt(3*A)))*2*(2*h/A)**(I(1)/6)+iv.sqrt(2*iv.pi*h/B)
    p=2*iv.exp(-delta**2/(T*h)); pW=p*W
    arc=(iv.pi/4)*iv.exp(abs(a)/(2*h)-(T/h)*(I(19)/(64*448)-delta**4/7))
    # Bound difference of the Q=1 and Q=1/2 reduced integrals using unshifted damping.
    c=T*3/128
    assert abs(a).b<=c.a
    ext=iv.exp(abs(a)/(8*h)-c/(16*h))*(1+(3*T*M+T*d*g/h)/2)
    # endpoint tilted law: actual spectral mass and diagonal bound.
    f=2*K*d*delta**2*g; assert f.b<1
    sig=[iv.sqrt(h*x/(1-f)) for x in G]
    D=iv.exp(K*d*delta**2*sum(G)/(1-f))
    # Build polynomial in scaled t=s*x and common scalar standard normal |Z|.
    X=[]
    for j,cj in enumerate([1,3,8,22],1):
        X.append({(4-j-l,j+l):cj*d/h*math.comb(4-j,l)*s**(4-j-l)*sum(vv**(4-j-l)*ss**(j+l) for vv,ss in zip(v,sig)) for l in range(5-j)})
    x1,x2,x3,x4=X; x11=polymul(x1,x1)
    P=x4
    for p0 in [scale(polymul(x1,x3),4),scale(polymul(x2,x2),3),scale(polymul(x11,x2),6),polymul(x11,x11)]: P=polyadd(P,p0)
    # Moment bounds for each product: generalized Holder over all Gaussian factors.
    def expect(p0):
        out={}
        for (i,j),v0 in p0.items(): out[i]=out.get(i,I(0))+v0*gm(j)
        return out
    P=expect(P)
    # Upper Riemann sums on scaled [0,1/s] (round upward by one panel).
    # H is the signed or unsigned quadratic envelope with quartic damping 3*T*c4/8.
    # Integrating the positive envelope slightly beyond the true endpoint is safe.
    mesh=I(1)/100; steps=math.ceil(100*upper(1/s))+1
    assert (steps*mesh).a >= (1/s).b
    b=3*T*c4/8
    signed=-a/(2*iv.sqrt(2*h))
    def envelope(lo,hi):
        if variant=="baseline": return iv.exp(k*hi**2/2-b*lo**4)
        return iv.exp(max((signed*lo**2).b,(signed*hi**2).b)-b*lo**4)
    E=I(0); E2=I(0)
    for j in range(steps):
        lo=j*mesh; hi=(j+1)*mesh; H=envelope(lo,hi)
        E+=mesh*H*sum(v0*hi**i for i,v0 in P.items())
    # For omitted polynomial use Cauchy-Schwarz (square polynomial then moment).
    P2sq=expect(polymul(polyadd({(0,0):I(1)},scale(polyadd(x2,x11),I(1)/2)),polyadd({(0,0):I(1)},scale(polyadd(x2,x11),I(1)/2))))
    for j in range(steps):
        lo=j*mesh; hi=(j+1)*mesh; H=envelope(lo,hi)
        E2+=mesh*H*iv.sqrt(sum(v0*hi**i for i,v0 in P2sq.items()))
    B4=D*2*E/24; Bpoly=2*iv.sqrt(p)*E2
    items={'k':k,'T':T,'tilt':f,'D':D,'scaled_real_tail':tail/s,'scaled_discard':pW/s,'scaled_arc':arc/s,'scaled_extension':ext/s,'scaled_B4':B4,'scaled_Bpoly':Bpoly,'scaled_total':(tail+pW+arc+ext)/s+B4+Bpoly}
    print('N,h,eta,variant',N,hs,es,variant,flush=True)
    for name,value in items.items(): print(name,value,flush=True)
    return items

if __name__=='__main__':
    print(sys.version, mp.__version__,flush=True)
    for variant in ['baseline','sharp']:
        for N,hs,es in [(4,'1e-5','0'),(80,'1e-5','0.001'),(80,'1e-5','-0.001'),(80,'1e-6','0.0003')]: run(N,hs,es,variant)
