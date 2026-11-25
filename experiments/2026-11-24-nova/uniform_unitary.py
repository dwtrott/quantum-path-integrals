"""All-mesh envelope for the unitary bound. Uses mpmath intervals throughout.
Run from this directory: OPENBLAS_NUM_THREADS=1 python uniform_unitary.py
"""
from mpmath import iv
from green_a1 import linear

iv.dps=40
I=iv.mpf

def hull(x,y): return I([min(x.a,y.a),max(x.b,y.b)])
def mag(x): return max(abs(x.a),abs(x.b))
def pad(x,b): return x+I([-b.b,b.b])
def force(x,jac=False):
    def at(v): return 3*(1-v*v if jac else v)*iv.exp(-v*v/2)
    out=hull(at(x.a),at(x.b))
    if jac:
        k=iv.sqrt(I(3))
        for v in (-k,k):
            if x.a<=v.b and v.a<=x.b: out=hull(out,-6*iv.exp(-I(3)/2))
        if x.a<=0<=x.b: out=hull(out,I(3))
    else:
        for v in (-1,1):
            if x.a<=v<=x.b: out=hull(out,at(I(v)))
    return out

def step(q,v,h,d,c=None):
    # d is an upper bound on the unknown mesh spacing, h is the tube slab length.
    def rhs(Q,V):
        b=((mag(V)+I(3)*d*max(I(1),mag(Q))/2)/(1-I(3)*d*d/2)).b
        nodes=pad(Q,d*b)
        f=force(nodes) if c is None else c*nodes
        return pad(V,d*mag(f)/2),f
    Q=pad(q,2*h*(1+mag(v)+mag(q))); V=pad(v,8*h*(1+mag(q)+mag(v)))
    b,f=rhs(Q,V); unit=I([0,h.b])
    Q1,V1=q+unit*b,v+unit*f
    assert Q.a<Q1.a and Q1.b<Q.b and V.a<V1.a and V1.b<V.b
    b,f=rhs(Q1,V1)
    return q+h*b,v+h*f,Q1

def large(n0=128,steps=1024,blocks=16):
    d=I(2)/n0; h=I(2)/steps
    x,v=I(-2),I(2); xs=[]
    for j in range(steps):
        x,v,tube=step(x,v,h,d); xs.append(tube)
    reach=(steps+n0-1)//n0+1
    cs=[]
    for j in range(steps):
        boxes=xs[max(0,j-reach):min(steps,j+reach+1)]
        cs.append(force(I([min(x.a for x in boxes),max(x.b for x in boxes)]),True))
    def shoot(cs):
        a,b,av,bv=I(1),I(0),I(0),I(1); out=[]
        for c in cs:
            a,av,A=step(a,av,h,d,c); b,bv,B=step(b,bv,h,d,c)
            out.append((A**2+B**2).b)
        return out,((av-b)**2+(bv+a)**2).b
    l,w2=shoot(cs); u,_=shoot(cs[::-1]); u.reverse()
    assert steps%blocks==0
    width=steps//blocks; mass=I(2)/blocks+d
    S4=S6=prefix=I(0)
    for k in range(blocks):
        L=max(l[k*width:(k+1)*width]); U=max(u[k*width:(k+1)*width])
        S4+=mass*L*U
        f,b=mass*L*iv.sqrt(L),mass*U*iv.sqrt(U)
        S6+=b*(f+2*prefix); prefix+=f
    C=(9*S4/32+15*23*S6/576)*iv.sqrt(iv.sqrt(w2))/iv.sqrt(2)
    print('large',n0,steps,blocks,'W2 upper',w2,'S4 upper',S4.b,'S6 upper',S6.b,'C upper',C.b,flush=True)
    return C.b

if __name__=='__main__':
    finite=[linear(n,True)[4].b for n in range(1,128)]
    print('finite N<=15 max C upper',max(finite[:15]),'N=',1+finite.index(max(finite[:15])),flush=True)
    print('finite 16<=N<=127 max C upper',max(finite[15:]),'N=',1+finite.index(max(finite[15:])),flush=True)
    C=large()
    assert max(finite)<560 and max(finite[15:])<68 and C<68
