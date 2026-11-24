"""Interval check for correspondence 049. Standard-library Python, IEEE binary64."""
from decimal import Decimal, localcontext
from math import inf, nextafter
from functools import lru_cache
import sys
def down(x): return nextafter(x,-inf)
def up(x): return nextafter(x,inf)
class I:
    def __init__(self,l,h=None):
        self.l=l
        self.h=l if h is None else h
    def __add__(a,b):
        if not isinstance(b,I): b=I(b)
        return I(down(a.l+b.l),up(a.h+b.h))
    __radd__=__add__
    def __neg__(a): return I(-a.h,-a.l)
    def __sub__(a,b): return a+(-b)
    def __rsub__(a,b): return -a+b
    def __mul__(a,b):
        if not isinstance(b,I): b=I(b)
        p=[a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h]
        return I(down(min(p)),up(max(p)))
    __rmul__=__mul__
    def __truediv__(a,b):
        if not isinstance(b,I): b=I(b)
        assert b.l>0
        return a*I(down(1/b.h),up(1/b.l))
    def mag(a): return max(abs(a.l),abs(a.h))
    def hull(a,b): return I(min(a.l,b.l),max(a.h,b.h))
    def sq(a): return I(max(0,down(min(abs(a.l),abs(a.h))*min(abs(a.l),abs(a.h)))) if (a.l>0 or a.h<0) else 0,up(a.mag()*a.mag()))
    def pad(a,e): return a+I(-e,e)
@lru_cache(maxsize=100000)
def exp_point(x):
    with localcontext() as ctx:
        ctx.prec=48
        v=Decimal.from_float(x).exp()
        return I(down(float(v.next_minus())),up(float(v.next_plus())))
def sqrt_interval(x):
    assert x.l >= 0
    with localcontext() as ctx:
        ctx.prec = 48
        lower = Decimal.from_float(x.l).sqrt().next_minus()
        upper = Decimal.from_float(x.h).sqrt().next_plus()
        return I(down(float(lower)), up(float(upper)))

def exp(x): return I(exp_point(x.l).l,exp_point(x.h).h)
def endpoint(x,jac):
    z=I(x)
    w=exp(-.5*z.sq())
    return 3*(1-z.sq() if jac else z)*w
def force(x,jac=False):
    v=endpoint(x.l,jac).hull(endpoint(x.h,jac))
    if jac:
        for lo,hi in [(-1.733,-1.732),(1.732,1.733)]:
            if x.l<=hi and x.h>=lo: v=v.hull(-6*exp_point(-1.5))
        if x.l<=0<=x.h: v=v.hull(I(3))
    else:
        for c in [-1,1]:
            if x.l<=c<=x.h: v=v.hull(endpoint(c,False))
    return v
def tube_step(q,v,h,d,curvature=None):
    # Linear interpolants of the Verlet nodal positions and velocities.
    def rhs(Q,V):
        speed=((I(V.mag())+1.5*d*max(1,Q.mag()))/(1-1.5* d.sq()).l).h
        nodes=Q.pad((d*speed).h)
        f=force(nodes) if curvature is None else curvature*nodes
        return V.pad((.5*d*f.mag()).h),f
    tubeq=q.pad((2*h*(1+v.mag()+q.mag())).h)
    tubev=v.pad((8*h*(1+q.mag()+v.mag())).h)
    b,c=rhs(tubeq,tubev)
    unit=I(0,h.h)
    def inside(x,y): return y.l<x.l and x.h<y.h
    assert inside(q+unit*b,tubeq) and inside(v+unit*c,tubev)
    # Refine the tube once before propagating.
    Q,V=q+unit*b,v+unit*c
    b,c=rhs(Q,V)
    return q+h*b,v+h*c,Q
def sum_up(values):
    v=I(0)
    for x in values: v+=x
    return v
def large(n0,steps):
    d=I(2)/n0
    h=I(2)/steps
    x,v=I(-2),I(2)
    xs=[]
    for j in range(steps):
        x,v,tube=tube_step(x,v,h,d)
        xs.append(tube)
    # Curvature at both ends of every discrete segment crossing each slab.
    reach=(d/h.l).h
    radius=int(reach)+1
    cs=[]
    for j in range(steps):
        boxes=xs[max(0,j-radius):min(steps,j+radius+1)]
        cs.append(force(I(min(x.l for x in boxes),max(x.h for x in boxes)),True))
    def solve(coeff):
        a,b,av,bv=I(1),I(0),I(0),I(1)
        vals=[]
        for c in coeff:
            a,av,aa=tube_step(a,av,h,d,c)
            b,bv,bb=tube_step(b,bv,h,d,c)
            vals.append((aa.sq()+bb.sq()).h)
        return vals, (av-b).sq()+(bv+a).sq()
    left,den=solve(cs)
    right,_=solve(reversed(cs))
    # Pointwise Green bound (i<=j); tighter than multiplying global maxima.
    running=0.
    bound=0.
    for l,r in zip(left,reversed(right)):
        running=max(l,running)
        bound=max(bound,(I(running)*r/den).h)
    # d * (N+1) = 2+d
    hs2=((2+d).sq()*bound).h
    blocks=32
    assert steps%blocks==0
    width=steps//blocks
    lv=[max(left[i:i+width]) for i in range(0,steps,width)]
    rr=right[::-1]
    rv=[max(rr[i:i+width]) for i in range(0,steps,width)]
    weight=(I(2)/blocks+d).sq()
    total=I(0)
    for j in range(blocks):
        total+=weight*rv[j]*(lv[j]+2*sum_up(lv[:j]))/den
    hs2=min(hs2,total.h)
    print("large",n0,steps,"Wronskian squared lower",den.l,"d^2 HS^2 upper",hs2,flush=True)
    return hs2
def verlet(q,v,d,f0,at_next):
    mid=v+.5*d*f0
    q=q+d*mid
    return q,mid+.5*d*at_next(q)
def finite(n):
    d=I(2)/n
    x,v=I(-2),I(2)
    c=[force(x,True)]
    for j in range(n):
        x,v=verlet(x,v,d,force(x),force)
        c.append(force(x,True))
    def solve(cs):
        a,b,av,bv=I(1),I(0),I(0),I(1)
        out=[I(1)]
        for c0,c1 in zip(cs,cs[1:]):
            a,av=verlet(a,av,d,c0*a,lambda x:c1*x)
            b,bv=verlet(b,bv,d,c0*b,lambda x:c1*x)
            out.append(a.sq()+b.sq())
        return out,(av-b).sq()+(bv+a).sq()
    l,den=solve(c)
    r,_=solve(c[::-1])
    total,prefix=I(0),I(0)
    for lj,rj in zip(l,r[::-1]):
        total+=(lj+2*prefix)*rj
        prefix+=lj
    return (d.sq()*total/den).h
def strip():
    # |W'''|^2 depends only on x^2 and y^2. Enclose the entire strip.
    y2=I(0,(I(1)/5).sq().h)
    hi=0
    for j in range(6000):
        x2=(I(j,j+1)/1000).sq()
        poly=(x2+y2)*((3-x2+y2).sq()+4*x2*y2)
        hi=max(hi,(poly*exp(y2-x2)).h)
    # For x>=6, |z|(3+|z|^2)<=2*x^3; x^6 exp(-x^2) decreases.
    hi=max(hi,(4*I(6**6)*exp(y2-36)).h)
    assert hi<2.25
    return hi
if __name__=="__main__":
    assert sys.float_info.radix == 2 and sys.float_info.mant_dig == 53
    assert sys.float_info.rounds == 1
    assert (I(1.732).sq().h < 3 < I(1.733).sq().l)
    results=[finite(n) for n in range(1,128)]
    print("finite maximum d^2 HS^2",max(results),"N",1+results.index(max(results)))
    big=large(128,1024)
    assert max(max(results),big)<(I(19)/10).sq().l
    print("strip M3 squared upper",strip())
    D = I(19)/10
    pi_lower = (I(314)/100).l
    kappa = I(3)/10
    eta = kappa*16*D/I(pi_lower).sq()
    row = (4*sqrt_interval(D/2)
           + 6*sqrt_interval(I(2))*sqrt_interval(sqrt_interval(7*D)))/pi_lower
    assert eta.h < (I(93)/100).l and row.sq().h < 42
    # Convex chord: -log(1-x) <= 3x on [0,.93].
    assert (exp(3*I(93)/100)*7/100).l > 1
    fourth = exp(I(1)/50)*(3+I(8)/100+I(1)/625+6+I(6)/25+3)
    assert fourth.h < 14
    c1 = I(3)*14*2*42**2/8
    c2 = exp(I(38))*I(4.5).sq()*15/72*4*600**3
    assert (c1+c2).h < exp_point(61).l
    print("positivity product upper", eta.h, "row variance upper", row.sq().h)

