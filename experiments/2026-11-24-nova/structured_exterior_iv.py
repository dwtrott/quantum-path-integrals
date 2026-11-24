"""Conditional exterior bound for the structured contour; mpmath interval arithmetic.
OPENBLAS_NUM_THREADS=1 python structured_exterior_iv.py
"""
import numpy as np
from scipy.linalg import solve
from mpmath import iv, __version__ as mpversion
from structured_audit import setup
iv.dps=40

def I(x):return iv.mpf(str(x))
def Z(x):return iv.mpc(I(x.real),I(x.imag))
def upper(x):return float(x.b)
def certify(N):
    d,m,E,K,D,t,z,phi,grad=setup(N);n=N+1;dd=I(2)/N
    mm=[dd/2 if j in (0,N) else dd for j in range(n)]
    tt=[I('.02')+I('1.4')*(I(j)/N)*(1-I(j)/N) for j in range(n)]
    sigma=I('.8291');a=I('-1.6247');B=I('.4');eta=I('.3');rho=I('1e-8')
    # Exact real discrete classical target (interval enclosure).
    x=I(-2);p=I(2)+3*dd*x*iv.exp(-x*x/2)/2
    for j in range(1,n):
        x=x+dd*p;p=p+3*mm[j]*x*iv.exp(-x*x/2)
    qb=x-I('.6');pb=p+I('.3');zz=[Z(v) for v in z]
    def lap(v,j):return sum((v[j]-v[k])/dd for k in (j-1,j+1) if 0<=k<n)
    ff=[lap(zz,j)+3*mm[j]*zz[j]*iv.exp(-zz[j]**2/2) for j in range(n)]
    ff[0]+=1j*(zz[0]+2)+2;ff[-1]+=1j*(zz[-1]-qb)-pb
    hh=[I(1 if j in (0,N) else 2)/dd+3*mm[j]*(1-zz[j]**2)*iv.exp(-zz[j]**2/2)+(1j if j in (0,N) else 0) for j in range(n)]
    H=K+1j*D+np.diag(3*m*(1-z*z)*np.exp(-z*z/2));R=solve(H,np.eye(n,dtype=complex));RR=[[Z(x) for x in row] for row in R]
    # Ball in complex sup norm. R is a fixed approximate inverse, all entries interpreted as decimal rationals.
    resid=max(sum((RR[i][j]*ff[j] for j in range(n)),iv.mpc(0)).__abs__().b for i in range(n))
    defect=max(sum(abs((1 if i==j else 0)-RR[i][j]*hh[j]+sum(RR[i][k]/dd for k in (j-1,j+1) if 0<=k<n)) for j in range(n)).b for i in range(n))
    rn=max(sum(abs(v) for v in row).b for row in RR)
    variation=max((3*mm[j]*rho*(3*(abs(zz[j])+rho)+(abs(zz[j])+rho)**3)*iv.exp((abs(zz[j].imag)+rho)**2/2)).b for j in range(n))
    contraction=defect+rn*variation;assert contraction<1 and resid+contraction*rho<rho
    # R is nonsingular because ||I-R H(z0)||<1; fixed point gives actual stationary root.
    balls=[iv.mpc(v.real+iv.mpf([-1,1])*rho,v.imag+iv.mpf([-1,1])*rho) for v in zz]
    gq=[-3*mm[j]*v*iv.exp(-v**2/2) for j,v in enumerate(balls)]
    beta=B*(1+abs(a)/iv.sqrt(2));av=[(tt[j]+tt[j+1])/2 for j in range(N)]
    c=[B*abs(tt[j+1]-tt[j])/(sigma-max(tt[j],tt[j+1])) for j in range(N)]
    # Tridiagonal remaining form Ar = K T + T K + D(1-T²) - 2 eta E^T diag(av/d) E.
    off=[-2*(1-eta)*v/dd for v in av]
    diag=[2*tt[j]*(1 if j in (0,N) else 2)/dd+ (1-tt[j]**2 if j in (0,N) else 0)-2*eta*sum(av[k]/dd for k in (j-1,j) if 0<=k<N) for j in range(n)]
    pivot=diag[0]-I('.85')*mm[0];assert pivot>0;minpivot=pivot.a
    for j in range(1,n):
        pivot=diag[j]-I('.85')*mm[j]-off[j-1]**2/pivot;assert pivot>0;minpivot=min(minpivot,pivot.a)
    ell=[abs(gq[j].imag+tt[j]*gq[j].real)+(beta*tt[j] if j in (0,N) else 0) for j in range(n)]
    # SPD with nonpositive off-diagonals => inverse entrywise nonnegative. Verify a positive supersolution.
    mat=np.diag([upper(v) for v in diag])+np.diag([upper(v) for v in off],1)+np.diag([upper(v) for v in off],-1)
    s=solve(mat,np.array([upper(v)+.0001 for v in ell]));ss=[I(v) for v in s]
    assert all(v>0 for v in ss)
    for j in range(n):assert diag[j]*ss[j]+sum(off[min(j,k)]*ss[k] for k in (j-1,j+1) if 0<=k<n)>ell[j]
    L=iv.sqrt(sum(ell[j]*ss[j] for j in range(n)))
    E1=beta*sum(abs(v.real) for v in gq)+beta**2+3*sum(mm[j]*(iv.exp((abs(v.imag-tt[j]*v.real)+beta)**2/(2*(1-tt[j]**2)))+abs(iv.exp(-v**2/2))) for j,v in enumerate(balls))+sum(c[j]**2/(4*eta*av[j]*dd) for j in range(N))
    # sufficient exterior radius 8, with margin uniform over the certified stationary ball.
    assert L<I('.75') and E1<I('19.5'); lower=I('.5')*I('.85')*64-I('.75')*iv.sqrt(I('.85'))*8-I('19.5');assert lower>0
    print('N',N,'mpmath',mpversion,'resid',resid,'contraction',contraction,'minpivot',minpivot,'L',L,'E',E1,'h at exterior lower',lower,flush=True)

if __name__=='__main__':
    for N in [16,32,64,128]:certify(N)
