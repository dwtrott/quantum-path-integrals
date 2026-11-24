"""Two-sided Hermite majorant for 0<h<=1/20, including off-grid h."""
from math import factorial,prod
from uniform_049 import I,sqrt_interval,force,verlet,tube_step
def budgets(exposures):
    # Terms of the positive exponential generating function at orders <=4.
    deriv=[0,I(4)/5,I(3)/8,I(4)/25,I(1)/16]
    out=[]
    for degree in range(5):
        total=I(0)
        for p in range(5):
            for q in range(3):
                for r in range(2):
                    for s in range(2):
                        ns=(p,q,r,s)
                        if p+2*q+3*r+4*s!=degree:continue
                        D=degree+2*sum(ns)
                        term=sqrt_interval(I(prod(range(1,2*D,2)))/(2**D))
                        for j,n in enumerate(ns,1):
                            for k in range(1,n+1):term=term*deriv[j]*exposures[j-1]/k
                        total+=term
        out.append(total)
    return out
def combine(l,r,norm):
    L,R=budgets(l),budgets(r)
    h=I(1)/20
    return (norm*(L[4]+R[4]+h.sq()*L[4]*R[4]+L[1]*R[3]+L[2]*R[2]+L[3]*R[1]+h*L[3]*R[3])).h
def evaluate(left,right,mleft,mright,den):
    def exposure(rows,mass):
        ans=[I(0) for _ in range(4)]
        for row,t in zip(rows,mass):
            b=sqrt_interval(I(row))
            v=b*b*b*t
            for j in range(4):
                ans[j]+=v
                v*=b
        return ans
    return combine(exposure(left,mleft),exposure(right,mright),sqrt_interval(sqrt_interval(den)/2))
def fixed(n):
    d=I(2)/n;x,v=I(-2),I(2);cs=[force(x,True)]
    for _ in range(n):
        x,v=verlet(x,v,d,force(x),force);cs.append(force(x,True))
    def shoot(cs):
        a,b,av,bv=I(1),I(0),I(0),I(1);res=[1]
        for c,e in zip(cs,cs[1:]):
            a,av=verlet(a,av,d,c*a,lambda x:e*x)
            b,bv=verlet(b,bv,d,c*b,lambda x:e*x)
            res.append((a.sq()+b.sq()).h)
        return res,(av-b).sq()+(bv+a).sq()
    l,den=shoot(cs);r,_=shoot(cs[::-1])
    # Split the kick clock at a Voronoi boundary, or inside a kick.
    cut=I(1) if n in (1,3) else I(3)/2 if n==2 else I(5)/4
    ml=[];mr=[]
    for j in range(n+1):
        lo=I(0) if j==0 else d*(j-.5)
        hi=I(2) if j==n else d*(j+.5)
        ml.append(I(max(0,(I(min(hi.h,cut.h))-lo).h)))
        mr.append(I(max(0,(hi-I(max(lo.l,cut.l))).h)))
    return evaluate(l,r[::-1],ml,mr,den)
def large(n0=128,steps=1024):
    D=I(2)/n0;h=I(2)/steps
    x,v=I(-2),I(2);xs=[]
    for _ in range(steps):
        x,v,box=tube_step(x,v,h,D);xs.append(box)
    reach=(steps+n0-1)//n0+1
    def near(seq,j):return seq[max(0,j-reach):j+reach+1]
    cs=[force(I(min(z.l for z in near(xs,j)),max(z.h for z in near(xs,j))),True) for j in range(steps)]
    def shoot(cs):
        a,b,av,bv=I(1),I(0),I(0),I(1);res=[]
        for c in cs:
            a,av,A=tube_step(a,av,h,D,c);b,bv,B=tube_step(b,bv,h,D,c)
            res.append((A.sq()+B.sq()).h)
        return res,(av-b).sq()+(bv+a).sq()
    l,den=shoot(cs);r,_=shoot(cs[::-1]);r.reverse()
    ml=[h if 8*j<5*steps else I(0) for j in range(steps)]
    mr=[h if 8*j>=5*steps else I(0) for j in range(steps)]
    return evaluate([max(near(l,j)) for j in range(steps)],[max(near(r,j)) for j in range(steps)],ml,mr,den)
if __name__=='__main__':
    v=[fixed(n) for n in range(1,128)]
    fine=large()
    print('finite max',max(v),'N',1+v.index(max(v)),'N>=16',max(v[15:]),'large',fine)
    assert max(max(v),fine)<20068 and max(max(v[15:]),fine)<8707
