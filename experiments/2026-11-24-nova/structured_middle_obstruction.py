"""Interval local bound and exact Lagrangian-gap witnesses for correspondence 039.
OPENBLAS_NUM_THREADS=1 python structured_middle_obstruction.py
"""
import numpy as np
import scipy
import mpmath
from mpmath import iv
from scipy.linalg import solve
from structured_audit import setup
iv.dps=160

def I(x): return iv.mpf(str(x))
def Z(x): return iv.mpc(I(x.real),I(x.imag))
def root_ball(N,z,K,D,m,d,rho):
    n=N+1; mm=np.array([float(v.a) for v in m])
    H=K+1j*D+np.diag(3*mm*(1-z*z)*np.exp(-z*z/2))
    R=solve(H,np.eye(n,dtype=complex));RR=[[Z(v) for v in row] for row in R]
    zz=[Z(v) for v in z]
    x=I(-2);p=I(2)+3*d*x*iv.exp(-x*x/2)/2
    for j in range(1,n):
        x=x+d*p;p=p+3*m[j]*x*iv.exp(-x*x/2)
    qb=x-I('.6');pb=p+I('.3')
    ff=[sum((zz[j]-zz[k])/d for k in (j-1,j+1) if 0<=k<n)+3*m[j]*zz[j]*iv.exp(-zz[j]**2/2) for j in range(n)]
    ff[0]+=1j*(zz[0]+2)+2;ff[-1]+=1j*(zz[-1]-qb)-pb
    hh=[I(1 if j in (0,N) else 2)/d+3*m[j]*(1-zz[j]**2)*iv.exp(-zz[j]**2/2)+(1j if j in (0,N) else 0) for j in range(n)]
    residual=max(abs(sum(RR[i][j]*ff[j] for j in range(n))).b for i in range(n))
    defect=max(sum(abs((1 if i==j else 0)-RR[i][j]*hh[j]+sum(RR[i][k]/d for k in (j-1,j+1) if 0<=k<n)) for j in range(n)).b for i in range(n))
    rn=max(sum(abs(v) for v in row).b for row in RR)
    variation=max((3*m[j]*rho*(3*(abs(zz[j])+rho)+(abs(zz[j])+rho)**3)*iv.exp((abs(zz[j].imag)+rho)**2/2)).b for j in range(n))
    contraction=defect+rn*variation;assert contraction<1 and residual+contraction*rho<rho
    print('root',N,'residual upper',float(residual.b),'contraction upper',float(contraction.b),flush=True)
    return zz

def certify(N):
    _,_,_,K,D,_,z,_,_=setup(N)
    d=I(2)/N;n=N+1;rho=I('1e-8');s=I('.8291');a=I('-1.6247');B=I('.4')
    m=[d/2 if j in (0,N) else d for j in range(n)]
    zc=root_ball(N,z,K,D,m,d,rho)
    zz=[v+iv.mpc(iv.mpf([-1,1])*rho,iv.mpf([-1,1])*rho) for v in zc]
    # Center Hessian. Root uncertainty charged separately in Taylor bound.
    h=[3*m[j]*(1-zc[j]**2)*iv.exp(-zc[j]**2/2)+(1j if j in (0,N) else 0) for j in range(n)]
    u=1+1j*s;b=1j*a/2;total=sum(h);lam=I('.65')
    # C=u I+b 1 m^T, K 1=0. Lower triangle of Im(C^T H C)-lambda M.
    A=[]
    for i in range(n):
        row=[]
        for j in range(i+1):
            val=(u*b*(h[i]*m[j]+h[j]*m[i])+b*b*total*m[i]*m[j]).imag
            if i==j:val+=(u*u*(h[i]+I(1 if i in (0,N) else 2)/d)).imag-lam*m[i]
            elif i==j+1:val+=(-u*u/d).imag
            row.append(val)
        A.append(row)
    if N<=64:
        L=[[I(0)]*n for _ in range(n)];piv=[]
        for i in range(n):
            for j in range(i):L[i][j]=(A[i][j]-sum(L[i][k]*piv[k]*L[j][k] for k in range(j)))/piv[j]
            p=A[i][i]-sum(L[i][k]**2*piv[k] for k in range(i));assert p>0,(N,i,p);piv.append(p)
        norm2=1+max((s*s).b,((s+a)**2).b)
        for rr in ['.02','.03','.04','.05','.06']:
            r=I(rr);assert s*r<B and iv.sqrt(2)*r<B
            ymax=(s+abs(a))*r;w=iv.mpc(iv.mpf([-1,1])*r,iv.mpf([-1,1])*ymax)
            rad=iv.sqrt(r*r+ymax*ymax)
            third=max((3*abs((3*(v+w)-(v+w)**3)*iv.exp(-(v+w)**2/2))).b for v in zz)
            margin=lam-norm2*third*(rho+rad/3)
            print('local',N,'r',rr,'2h/||v||_M^2 lower',float(margin.a),'pivot lower',float(min(p.a for p in piv)),flush=True)
            if rr=='.04':assert margin>I('.19')
    # Exact conditional objective at q=0 but v_j= +/- .3 (infeasible mean).
    # All clips inactive, kinetic differences vanish. Stationarity subtracts linear part.
    vals=[]
    for sign in [1,-1]:
        v=sign*I('.3');assert abs(s*v)<B
        w=u*v
        W=lambda x:iv.exp(-x*x/2)
        value=(1-s*s)*v*v+3*sum(m[j]*(W(zz[j])-W(zz[j]+w)-zz[j]*W(zz[j])*w).imag for j in range(n))
        assert value<I('-.06');vals.append(value)
        print('dual witness',N,sign,'lower',float(value.a),'upper',float(value.b),flush=True)
    average=(vals[0]+vals[1])/2;assert average<I('-.11')
    print('dual supremum upper',N,float(average.b),flush=True)
if __name__=='__main__':
    print('versions',np.__version__,scipy.__version__,mpmath.__version__)
    for N in [16,32,64,128]:certify(N)
