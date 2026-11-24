"""Interval Taylor certificate with a rank-two inertia reduction.
Run OPENBLAS_NUM_THREADS=1 python structured_local_rank.py.
"""
import os,sys,time
import numpy as np
import scipy,mpmath
from mpmath import iv
from structured_audit import setup
from structured_middle_obstruction import root_ball,I
iv.dps=45

def max_derivative(z,r,yr,rho,panels):
    # z has point-center uncertainty from decimal representation only.
    # Rectangle includes every z_actual + tau*w for tau in [0,tau_max].
    xr=r+rho; yi=yr+rho
    hi=I(0)
    for k in range(panels):
        x=z.real+xr*iv.mpf([I(2*k)/panels-1,I(2*(k+1))/panels-1])
        for l in range(panels):
            y=z.imag+yi*iv.mpf([I(2*l)/panels-1,I(2*(l+1))/panels-1])
            # |3(z^3-3z)e^(-z^2/2)| = 3|z| |z^2-3| exp((y^2-x^2)/2).
            x2=x**2;y2=y**2
            mag=3*iv.sqrt((x2+y2)*((x2-y2-3)**2+4*x2*y2))*iv.exp((y2-x2)/2)
            hi=max(hi,mag.b)
    return hi

def bound_errors(zz,r,s,a,rho,panels=6,radial=8):
    mean=I(os.environ["MEAN_RADIUS"]) if "MEAN_RADIUS" in os.environ else r
    assert mean>=0 and mean<=r
    yr=s*r+abs(a)*mean; L=iv.sqrt(r*r+yr*yr)
    out=[]
    for z in zz:
        e=rho*max_derivative(z,I(0),I(0),rho,1)
        for k in range(radial):
            l=I(k)/radial;h=I(k+1)/radial
            weight=((1-l)**3-(1-h)**3)/3
            e+=L*weight*max_derivative(z,h*r,h*yr,rho,panels)
        out.append(e.b)
    return out

def inertia_certificate(h,m,err,d,s,a,lam):
    n=len(m);u=1+1j*s;b=1j*a/2
    # B-lambda M = T + U J U^T, U=[m,g], J=[[beta,1],[1,0]].
    diag=[2*s*(1 if j in (0,n-1) else 2)/d+(u*u*h[j]).imag-(1+s*s)*m[j]*err[j]-lam*m[j] for j in range(n)]
    g=[(u*b*h[j]).imag-s*a/2*m[j]*err[j] for j in range(n)]
    beta=(b*b*sum(h)).imag-a*a/4*sum(m[j]*err[j] for j in range(n))
    off=-2*s/d; p=[diag[0]]
    if not (p[0]>0 or p[0]<0):return None
    for j in range(1,n):
        v=diag[j]-off**2/p[-1]
        if not (v>0 or v<0):return None
        p.append(v)
    neg=sum(v<0 for v in p)
    def forward(v):
        f=[v[0]]
        for j in range(1,n):f.append(v[j]-off/p[j-1]*f[-1])
        return f
    fm=forward(m);fg=forward(g)
    # -J^{-1} = [[0,-1],[-1,beta]].  Use T=LDL^T to get U^T T^-1 U.
    A=-sum(fm[j]**2/p[j] for j in range(n));B=-1-sum(fm[j]*fg[j]/p[j] for j in range(n));C=beta-sum(fg[j]**2/p[j] for j in range(n))
    det=A*C-B**2
    # For this problem T has one negative eigenvalue. The 2x2 Schur block must be positive definite.
    ok=neg==1 and A>0 and det>0
    return ok,neg,float(A.a),float(det.a),float(min(abs(v).a for v in p))

def certify(N,r,panels=6,radial=8):
    start=time.perf_counter();d0,m0,E,K,D,t,z,*_=setup(N);d=I(2)/N;rho=I('1e-8');s=I('.8291');a=I('-1.6247');B=I('.4');r=I(r);n=N+1
    assert s*r<B and iv.sqrt(2)*r<B
    m=[d/2 if j in (0,N) else d for j in range(n)]
    zz=root_ball(N,z,K,D,m,d,rho)
    h=[3*m[j]*(1-zz[j]**2)*iv.exp(-zz[j]**2/2)+(1j if j in (0,N) else 0) for j in range(n)]
    err=bound_errors(zz,r,s,a,rho,panels,radial)
    results=[]
    for lam in ['0','.025','.05','.1','.15','.2']:
        v=inertia_certificate(h,m,err,d,s,a,I(lam));results.append((lam,v))
    print('N',N,'r',r,'mean_radius',os.environ.get('MEAN_RADIUS','r'),'panels',panels,'radial',radial,'max err',float(max(err)),'seconds',time.perf_counter()-start,'results',results,flush=True)
    return results

if __name__=='__main__':
    print('versions',np.__version__,scipy.__version__,mpmath.__version__)
    r=sys.argv[1] if len(sys.argv)>1 else '.12'; panels=int(sys.argv[2]) if len(sys.argv)>2 else 6; radial=int(sys.argv[3]) if len(sys.argv)>3 else 8
    for N in map(int,sys.argv[4:] or [16,32,64,128]):certify(N,r,panels,radial)
