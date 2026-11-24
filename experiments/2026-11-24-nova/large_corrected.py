"""All large discrete meshes, via a second-order interpolation of Verlet nodes."""
import math
from uniform_049 import I, exp, force, tube_step
from classical_scalar import Z

def hull(items): return I(min(s.l for s in items),max(s.h for s in items))
def deriv(x):
    w=3*exp(-.5*x.sq())
    return force(x),force(x,True),w*x*(x.sq()-3),w*(-x.sq().sq()+6*x.sq()-3)

def all_large(start=256,steps=2048):
    D=I(2)/start; h=I(2)/steps
    reach=(steps+start-1)//start+2
    def coarse(seq,a,v,jac):
        out=[]
        for c in seq:
            old=v
            a,v,Q=tube_step(a,v,h,D,c if jac else None)
            # Same first-exit velocity box as in tube_step, enlarged using
            # |q_initial| <= mag(Q).
            V=old.pad((8*h*(1+Q.mag()+old.mag())).h)
            out.append((Q,V))
        return out
    xs=coarse(range(steps),I(-2),I(2),False)
    def near(seq,j):return seq[max(0,j-reach):min(steps,j+reach+1)]
    cur=[force(hull([q for q,v in near(xs,j)]),True) for j in range(steps)]
    coarse_sets=[xs,coarse(cur,I(1),I(0),True),coarse(cur,I(0),I(1),True)]
    q=[I(-2),I(1),I(0)]; v=[I(2),I(0),I(1)]
    fine=[]
    action=I(0)
    for j in range(steps):
        # Enclose both nodal endpoints and the entire polynomial segment.
        bounds=[near(seq,j) for seq in coarse_sets]
        pos=[hull([a for a,b in box]) for box in bounds]
        vel=[hull([b for a,b in box]) for box in bounds]
        f,c,cp,cpp=deriv(pos[0])
        f0=[f,c*pos[1],c*pos[2]]
        pos=[a.pad((D.sq()*F.mag()/8).h) for a,F in zip(pos,f0)]
        # q'=v_left+s F_left, with 0<=s<=D
        speed=[a+I(0,D.h)*F for a,F in zip(vel,f0)]
        f,c,cp,cpp=deriv(pos[0])
        G=[f,c*pos[1],c*pos[2]]
        G1=[c*speed[0]]+[cp*speed[0]*pos[k]+c*speed[k] for k in (1,2)]
        G2=[cp*speed[0].sq()+c*f0[0]]+[
            cpp*speed[0].sq()*pos[k]+cp*f0[0]*pos[k]
            +2*cp*speed[0]*speed[k]+c*f0[k] for k in (1,2)]
        e=[(D.sq()*g.mag()/2).h for g in G1]
        ev=[(D.sq()*g.mag()/4).h for g in G2]
        f,c,_,_=deriv(q[0]); initial=[f,c*q[1],c*q[2]]
        # Parametric interpolation shares q,v at the exact mesh nodes.
        dq=[(h*v[k]+.5*h.sq()*G[k].pad(ev[k])).pad((h*e[k]).h) for k in range(3)]
        dv=[(h*initial[k]+.5*h.sq()*G1[k]).pad((h*ev[k]).h) for k in range(3)]
        if j==0: assert speed[2].l>0
        def width(z): return I(z.h)-z.l
        dx=width(q[0])+h*width(v[0])+.5*h.sq()*width(G[0].pad(ev[0]))+2*h*e[0]+D.sq()*f0[0].mag()/8
        dp=width(v[0])+h*width(initial[0])+.5*h.sq()*width(G1[0])+2*h*ev[0]+D*f0[0].mag()/2+D.sq()*G1[0].mag()/2
        beta=I(speed[0].mag())+D*f0[0].mag()/2
        action+=h*((dp.sq()+3*dx.sq())/2+D.sq()*beta.sq()/4)
        fine.append([q[k]+I(0,h.h)*speed[k] for k in range(3)])
        q=[a+b for a,b in zip(q,dq)];v=[a+b for a,b in zip(v,dv)]
    y=Z(q[1],q[2]); vy=Z(v[1],v[2]); wr=vy+Z(0,1)*y
    assert all(t[2].l>0 for t in fine[1:]) and q[2].l>0
    # Phase orientation at the first slab follows b(0)=0, b'>0 there.
    lead=Z(1)/(y.root()*((Z(1)-Z(0,1)*vy/y)*.5).root())
    p1=p3=z4=z1=z3=Z()
    for j in range(steps):
        half=(steps+2*start-1)//(2*start)+1
        node=fine[max(0,j-half):min(steps,j+half+1)]
        x,a,b=[hull([v[k] for v in node]) for k in range(3)]
        l=Z(a,b)
        u=Z(v[2],q[2])*a-Z(v[1],q[1])*b
        V=3*exp(-.5*x.sq())
        v3=V*x*(3-x.sq());v4=V*(x.sq().sq()-6*x.sq()+3)
        g=l*u/wr
        w=g*v3
        z4+=g*g*(h*v4)
        z1+=(l*w*h+2*p1)*u*w*h
        z3+=(l.cube()*(h*v3)+2*p3)*u.cube()*(h*v3)
        p1+=l*w*h;p3+=l.cube()*(h*v3)
    a1=Z(0,1)*(z4/8+(9*z1/wr+6*z3/wr.cube())/72)
    return q[0],v[0],lead,a1,action.h
if __name__=="__main__":
    import sys
    q,p,L,a,action=all_large(*map(int,sys.argv[1:]))
    print(q.l,q.h,p.l,p.h,L.mid(),L.radius(),a.mid(),a.radius(),action)
