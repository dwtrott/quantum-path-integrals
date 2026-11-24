"""Validated scalar data for the real Verlet target (uniform_049 supplies intervals)."""
import math
from uniform_049 import I, exp, sqrt_interval

class Z:
    def __init__(self,r=0,i=0):
        self.r=r if isinstance(r,I) else I(r)
        self.i=i if isinstance(i,I) else I(i)
    def __add__(a,b):
        if not isinstance(b,Z): b=Z(b)
        return Z(a.r+b.r,a.i+b.i)
    __radd__=__add__
    def __neg__(a): return Z(-a.r,-a.i)
    def __sub__(a,b):return a+(-b)
    def __mul__(a,b):
        if not isinstance(b,Z): b=Z(b)
        return Z(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r)
    __rmul__=__mul__
    def conj(a):return Z(a.r,-a.i)
    def __truediv__(a,b):
        if not isinstance(b,Z): b=Z(b)
        n=b.r.sq()+b.i.sq()
        c=a*b.conj()
        return Z(c.r/n,c.i/n)
    def cube(a):return a*a*a
    def root(a):
        v=sqrt_interval((sqrt_interval(a.r.sq()+a.i.sq())+a.r)/2)
        return Z(v,a.i/(2*v))
    def mid(a):return complex((a.r.l+a.r.h)/2,(a.i.l+a.i.h)/2)
    def radius(a):return a.r.h-a.r.l+a.i.h-a.i.l+1e-14
def classical(n):
    d=I(2)/n
    x,p=I(-2),I(2)
    m=[d/2]+[d]*(n-1)+[d/2]
    left=[];cs=[]; v3=[];v4=[]
    y,v=Z(1),Z(0,1)
    action=I(-4)
    for j,t in enumerate(m):
        w=3*exp(-.5*x.sq())
        c=w*(1-x.sq());cs.append(c)
        v3.append(t*w*x*(3-x.sq()))
        v4.append(t*w*(x.sq().sq()-6*x.sq()+3))
        left.append(y)
        if j: assert y.i.l>0
        action-=t*w
        p+=t*x*w;v+=y*(t*c)
        if j<n:
            action+=.5*d*p.sq();x+=d*p;y+=v*d
    action-=p*x
    wr=v+Z(0,1)*y
    lead=Z(1)/(y.root()*((Z(1)-Z(0,1)*v/y)*.5).root())
    right=[];u,uv=Z(1),Z(0,1)
    for j in range(n,-1,-1):
        right.append(u);uv+=u*(m[j]*cs[j])
        if j:u+=uv*d
    right.reverse()
    p1=p3=z4=z1=z3=Z()
    for l,u,b,c in zip(left,right,v3,v4):
        g=l*u/wr
        z4+=g*g*c
        t=g*b
        z1+=(l*t+2*p1)*u*t
        z3+=(l.cube()*b+2*p3)*u.cube()*b
        p1+=l*t;p3+=l.cube()*b
    a1=Z(0,1)*(z4/8+(9*z1/wr+6*z3/wr.cube())/72)
    return x,p,action,lead,a1
