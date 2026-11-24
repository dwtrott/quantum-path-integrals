from poly_soft import *

def add(a,b):
    c=a.copy()
    for k,v in b.items(): c[k]=c.get(k,I(0))+v
    return c

def mul(a,b):
    c={}
    for (i,j),v in a.items():
        for (k,l),w in b.items(): c[i+k,j+l]=c.get((i+k,j+l),I(0))+v*w
    return c

def scale(a,c): return {k:v*c for k,v in a.items()}
def gm(k):
    if k%2==0: return I(math.prod(range(1,k,2)))
    return iv.sqrt(gm(k-1)*gm(k+1))

def budget(N,h,eta,meta,meshstr='0.01'):
    # Same sharp variant of 029, reuse the interval covariance diagonal from 028/026.
    T,a,g,G,v=meta; d=T/N; k=abs(a)/iv.sqrt(h); M=max(x.b for x in G)
    assert N>=4 and abs(eta).b<=I('0.5') and g.a>0 and min(x.a for x in G)>0
    delta=I(35)/1000; K=626; c4=I(1)/16; s=iv.sqrt(iv.sqrt(h)); A=T/512; B=T/128
    assert abs(a).b<=min(A.a,B.a)
    ct=A/(2*h); cg=B/(2*h); tail=2*iv.exp(-ct)/(6*ct)+2*iv.exp(-2*cg)/(2*cg*iv.sqrt(2))
    W=iv.exp(k**I('1.5')*h**I('-.25')/(3*iv.sqrt(3*A)))*2*(2*h/A)**(I(1)/6)+iv.sqrt(2*iv.pi*h/B)
    p=2*iv.exp(-delta**2/(T*h)); arc=(iv.pi/4)*iv.exp(abs(a)/(2*h)-(T/h)*(I(19)/(64*448)-delta**4/7))
    c=T*3/128; assert abs(a).b<=c.a
    ext=iv.exp(abs(a)/(8*h)-c/(16*h))*(1+(3*T*M+T*d*g/h)/2)
    f=2*K*d*delta**2*g; assert f.b<1
    sig=[iv.sqrt(h*x/(1-f)) for x in G]; D=iv.exp(K*d*delta**2*sum(G)/(1-f))
    X=[]
    for j,cj in enumerate([1,3,8,22],1):
        X.append({(4-j-l,j+l):cj*d/h*math.comb(4-j,l)*s**(4-j-l)*sum(vv**(4-j-l)*ss**(j+l) for vv,ss in zip(v,sig)) for l in range(5-j)})
    x1,x2,x3,x4=X; x11=mul(x1,x1); P=x4
    for p0 in [scale(mul(x1,x3),4),scale(mul(x2,x2),3),scale(mul(x11,x2),6),mul(x11,x11)]: P=add(P,p0)
    def expect(p0):
        out={}
        for (i,j),v0 in p0.items(): out[i]=out.get(i,I(0))+v0*gm(j)
        return out
    P=expect(P); P2=add({(0,0):I(1)},scale(add(x2,x11),I(1)/2)); P2sq=expect(mul(P2,P2))
    mesh=I(meshstr); steps=math.ceil(float((1/(s*mesh)).b))+1
    assert (steps*mesh).a >= (1/s).b
    b=3*T*c4/8; signed=-a/(2*iv.sqrt(2*h)); E=I(0); E2=I(0)
    # Evaluate both upper Riemann sums on the same mesh; Horner with nonnegative coefficients.
    def coeff(poly): return [poly.get(i,I(0)) for i in range(max(poly)+1)]
    def horner(C,x):
        y=C[-1]
        for c in reversed(C[:-1]): y=y*x+c
        return y
    C=coeff(P); C2=coeff(P2sq)
    for j in range(steps):
        lo=j*mesh; hi=(j+1)*mesh; H=iv.exp(max((signed*lo**2).b,(signed*hi**2).b)-b*lo**4)
        E+=mesh*H*horner(C,hi); E2+=mesh*H*iv.sqrt(horner(C2,hi))
    return (tail+p*W+arc+ext)/s+D*2*E/24+2*iv.sqrt(p)*E2

def full(N,hs,es,panels=24,meshstr='0.01'):
    t=time.perf_counter(); val,e,ep,tail,meta=run(N,hs,es,panels)
    B=budget(N,I(hs),I(es),meta,meshstr); R=B+e+ep+tail
    print('full',N,es,'mesh',meshstr,'stiff',B,'total disk before center rounding',R,'seconds',time.perf_counter()-t,flush=True)
    return R
if __name__=='__main__':
    for N,e in [(4,'0'),(80,'-0.001'),(80,'0.001'),(320,'0.001'),(1280,'0.001')]: full(N,'1e-5',e)
