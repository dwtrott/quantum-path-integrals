"""Validated complex shooting for Gaussian barrier coherent-state saddles.
Finite interval Taylor/Picard ODE enclosure, not a path-integral error bound.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from mpmath import iv
import mpmath,scipy,sys,time
iv.dps=40
I=iv.mpf; C=iv.mpc; J=C(0,1)

def encl(z):
    return C(I([z.real,z.real]),I([z.imag,z.imag]))
def box(z,r):
    rr=I([-r,r]); return z+C(rr,rr)
def inside(a,b):
    return a.real.a>b.real.a and a.real.b<b.real.b and a.imag.a>b.imag.a and a.imag.b<b.imag.b

def field(y):
    x,p,w,v=y; e=iv.exp(-x*x/2)
    return [p,3*x*e,v,3*(1-x*x)*e*w]

def coeff(y,n):
    """Ordinary Taylor coefficients of x,p,w,v and exp(-x²/2)."""
    x,p,w,v=([a] for a in y); e=[iv.exp(-x[0]*x[0]/2)]; a=[-x[0]*x[0]/2]; b=[]
    conv=lambda u,v,k:sum((u[j]*v[k-j] for j in range(k+1)),C(0))
    for k in range(n):
        if k:
            a.append(-conv(x,x,k)/2)
            e.append(sum((j*a[j]*e[k-j] for j in range(1,k+1)),C(0))/k)
        b.append((1 if k==0 else 0)+2*a[k])
        x.append(p[k]/(k+1));p.append(3*conv(x,e,k)/(k+1))
        w.append(v[k]/(k+1));v.append(3*sum((b[j]*conv(e,w,k-j) for j in range(k+1)),C(0))/(k+1))
    return [x,p,w,v]

def integrate(y,N=80,order=9):
    h=I(2)/N; dt=I([0,h.b]);
    for i in range(N):
        f=field(y); Y=[box(a, I('0.01').b+2*h.b*abs(b).b) for a,b in zip(y,f)]
        for k in range(20):
            Z=[a+dt*b for a,b in zip(y,field(Y))]
            if all(inside(a,b) for a,b in zip(Z,Y)):break
            Y=[box(a, I('0.001').b+abs(a-b).b) for a,b in zip(Y,Z)]
        else:raise RuntimeError(('Picard',i))
        # Taylor remainder coefficient at some point along enclosed trajectory;
        # complex rectangular integration bounds each real/imag remainder.
        A=coeff(y,order-1);B=coeff(Y,order)
        y=[sum((a[k]*h**k for k in range(order)),C(0))+b[order]*h**order for a,b in zip(A,B)]
    return y

def approx(z,p):
    def f(t,y):
        x,p=y;return [p,3*x*np.exp(-x*x/2)]
    a=solve_ivp(f,[0,2],[complex(z),complex(p)],rtol=2e-13,atol=2e-14,method='DOP853')
    assert a.success;return a.y[:,-1]

def verlet(z,p,w,v,N):
    d=I(2)/N
    def kick(z,p,w,v,h):
        e=iv.exp(-z*z/2)
        return p+h*3*z*e,v+h*3*(1-z*z)*e*w
    p,v=kick(z,p,w,v,d/2)
    for j in range(N):
        z=z+d*p;w=w+d*v
        p,v=kick(z,p,w,v,d/2 if j==N-1 else d)
    return z,p,w,v

def approx_verlet(z,p,N):
    d=2/N;p+=d/2*3*z*np.exp(-z*z/2)
    for j in range(N):
        z+=d*p;p+=(d/2 if j==N-1 else d)*3*z*np.exp(-z*z/2)
    return z,p

def certify_disk(z,rad,prop,Q,P):
    y=prop(z); F=y[1]-P+J*(y[0]-Q);D=y[3]+J*y[2]
    Z=box(z,rad.b);Y=prop(Z);Dbox=Y[3]+J*Y[2];a=C(D.real.mid,D.imag.mid)
    lhs=abs(F).b+rad*abs(Dbox-a).b;rhs=rad*abs(a).a
    assert lhs.b<rhs.a
    return F,Dbox,lhs,rhs

def run(N,order):
    start=time.perf_counter();q0=-2;p0=2
    qc,pc=approx(q0,p0)
    classical=integrate([C(q0),C(p0),C(0),C(0)],N,order)
    for dq,dp in [(0.4,0.4),(-0.6,0.3)]:
        Q=classical[0]+I(str(dq));P=classical[1]+I(str(dp))
        def residual(a):
            z=complex(*a);x,p=approx(z,p0+1j*(z-q0));r=p-(pc+dp)+1j*(x-(qc+dq));return [r.real,r.imag]
        sol=root(residual,[-2,0],tol=1e-11);assert np.linalg.norm(residual(sol.x))<1e-10
        z=encl(complex(*sol.x));rad=I('0.0000001')
        F,Dbox,lhs,rhs=certify_disk(z,rad,lambda z:integrate([z,C(p0)+J*(z-q0),C(1),J],N,order),Q,P)
        print('continuum target displacement',dq,dp,'center',sol.x.tolist(),'F',F,'Dbox',Dbox,'lhs',lhs,'rhs',rhs,'time',time.perf_counter()-start,flush=True)
        for n in (16,64,256):
            def residual_d(a):
                z=complex(*a);x,p=approx_verlet(z,p0+1j*(z-q0),n);r=p-(pc+dp)+1j*(x-(qc+dq));return [r.real,r.imag]
            sd=root(residual_d,sol.x,tol=1e-11);assert np.linalg.norm(residual_d(sd.x))<1e-10
            zd=encl(complex(*sd.x));F,Dbox,lhs,rhs=certify_disk(zd,rad,lambda z:verlet(z,C(p0)+J*(z-q0),C(1),J,n),Q,P)
            delta=abs(zd-z);scaled=delta*I(n*n)
            print('discrete',n,'fixed target',dq,dp,'center',sd.x.tolist(),'N^2 center displacement',float(scaled.mid),'lhs upper',float(lhs.b),'rhs lower',float(rhs.a),flush=True)
if __name__=='__main__':
    print('versions',np.__version__,scipy.__version__,mpmath.__version__,flush=True);run(int(sys.argv[1]) if len(sys.argv)>1 else 80,int(sys.argv[2]) if len(sys.argv)>2 else 9)
