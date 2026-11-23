"""Interval check of N=8, r=2 lattice using mpmath.iv; N=160 spectral bound.
Run: python pa_interval.py. mpmath interval backend is experimental; audit pending.
"""
import math
import mpmath as mp
iv=mp.iv
iv.dps=30
I=iv.mpf
J=iv.mpc(0,1)


def spectral(N,r):
    T=I(2); h=I(1)/10; tau=I(1)/N; t=(N-1)*tau; d=T/N
    lam=[(4*iv.sin(iv.pi*k/(2*N))**2-d*d)/d for k in range(1,N)]
    assert all(x.a>0 for x in lam)
    red=h*h*tau*tau*iv.exp(t)*sum(1/x**2 for x in lam[r:])/4
    # spectral bounds, no need to build dense covariance at large N
    B=h*2/N*sum(1/x for x in lam[:r])
    K=h*t*t*(I(2)/N)*sum(1/x for x in lam[r:])*B
    return lam,red,B,K


def quadrature_bound(r,B,K,t,step=I(1)/4,M=32,P=14):
    a=iv.pi**2/(step**2*(1+(P+2)*B)**2)
    v=iv.exp(-a); theta=1+2*v/(1-v)
    alias=iv.exp(t)*((1+K/2)*(theta**r-1)+K*iv.pi**2/step**2*r*theta**(r-1)*2*v*(1+v)/(1-v)**3)
    tail=t**(P+1)/math.factorial(P+1)/(1-t/(P+2))
    tay=tail*(1+K*r/2+(1+2*K/iv.exp(1))*(step+iv.sqrt(2*iv.pi))**r/iv.pi**(I(r)/2))
    trunc=iv.exp(t)*(1+2*K/iv.exp(1))*r*2*iv.exp(-step**2*M*M/2)/(step*M)*(step+iv.sqrt(2*iv.pi))**(r-1)/iv.pi**(I(r)/2)
    return alias+tay+trunc


def evaluate(N=8,r=2,M=32):
    assert r==2
    h=I(1)/10;tau=I(1)/N; t=(N-1)*tau;step=I(1)/4
    lam,red,B,K=spectral(N,r)
    V=[[iv.sqrt(I(2)/N)*iv.sin(iv.pi*j*k/N) for k in range(1,N)] for j in range(1,N)]
    H=[[sum(V[j][k]*V[l][k]/lam[k] for k in range(r,N-1)) for l in range(N-1)] for j in range(N-1)]
    D=[1+J*h*H[j][j] for j in range(N-1)]
    fac=[iv.exp(-iv.log(z)/2) for z in D]
    R=[[(1+J)*iv.sqrt(h/lam[k])*V[j][k] for k in range(r)] for j in range(N-1)]
    total=iv.mpc(0)
    for i in range(M+1):
        for l in range(M+1):
            x=i*step;y=l*step
            z=[row[0]*x+row[1]*y for row in R]
            W=[fac[j]*iv.exp(-z[j]*z[j]/(2*D[j])) for j in range(N-1)]
            wp=[-z[j]*W[j]/D[j] for j in range(N-1)]
            pair=sum(wp[j]*sum(H[j][k]*wp[k] for k in range(N-1)) for j in range(N-1))
            f=iv.exp(J*tau*sum(W))*(1+J*h*(J*tau)**2*pair/2)
            total+=(1 if i==0 else 2)*(1 if l==0 else 2)*iv.exp(-x*x-y*y)*f
    total*=step**r/iv.pi**(I(r)/2)
    quad=quadrature_bound(r,B.b,K.b,t,step,M)
    print('N',N,'lattice enclosure',total,'reduction bound',red,'quadrature bound',quad,'sum',red+quad)
    assert (red+quad).b<I('1e-5')


if __name__=='__main__':
    print('mpmath',mp.__version__)
    evaluate()
    N=160;r=2
    _,red,B,K=spectral(N,r)
    quad=quadrature_bound(r,B.b,K.b,I(N-1)/N)
    print('N',N,'reduction bound',red,'quadrature bound',quad,'sum',red+quad,'(lattice value not interval evaluated)')
    assert (red+quad).b<I('1e-5')
