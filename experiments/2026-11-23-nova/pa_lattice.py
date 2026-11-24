"""Nova: a finite lattice rule and analytic discretization bound for corrected PA.
Run: OPENBLAS_NUM_THREADS=1 python pa_lattice.py
Uses pa_benchmark.py. Floating point evaluation; roundoff not certified.
"""
import math
import numpy as np
from pa_benchmark import covariance, reduced


def lattice(V, lam, tau, h, r, step=.25, M=32, P=14):
    H = (V[:,r:]/lam[r:])@V[:,r:].T
    D = 1+1j*h*np.diag(H)
    b = np.sqrt(2*h)*(V[:,:r]/np.sqrt(lam[:r]))
    R = np.exp(1j*np.pi/4)*b
    B = np.max(np.sum(b*b,axis=1))/2
    t = len(lam)*abs(tau)
    K = h*t*t*np.max(abs(H))*B
    grid = step*np.arange(-M,M+1)
    assert r == 2, 'this implementation is for two retained modes'
    total = 0j
    for x in grid:
        q = np.column_stack((np.full(len(grid),x),grid))
        z = q@R.T
        W = np.exp(-z*z/(2*D))/np.sqrt(D)
        wp = -z/D*W
        f = np.exp(1j*tau*W.sum(axis=1))*(1+1j*h*(1j*tau)**2*np.sum((wp@H)*wp,axis=1)/2)
        total += np.sum(f*np.exp(-np.sum(q*q,axis=1)))
    total *= step**r/np.pi**(r/2)
    a = np.pi**2/(step**2*(1+(P+2)*B)**2)
    v = np.exp(-a)
    theta = 1+2*v/(1-v)
    alias = np.exp(t)*((1+K/2)*np.expm1(r*np.log1p(2*v/(1-v)))+K*np.pi**2/step**2*r*theta**(r-1)*2*v*(1+v)/(1-v)**3)
    tail = t**(P+1)/math.factorial(P+1)/(1-t/(P+2))
    tay = tail*(1+K*r/2+(1+2*K/np.e)*(step+np.sqrt(2*np.pi))**r/np.pi**(r/2))
    trunc = np.exp(t)*(1+2*K/np.e)*r*2*np.exp(-step**2*M*M/2)/(step*M)*(step+np.sqrt(2*np.pi))**(r-1)/np.pi**(r/2)
    return total, dict(B=B,K=K,alias=alias,taylor=tay,trunc=trunc,bound=alias+tay+trunc)


if __name__=='__main__':
    for N in (8,160):
        V,lam,_=covariance(N,2.); tau=1/N
        val,err=lattice(V,lam,tau,.1,2)
        gh,red=reduced(V,lam,tau,.1,2,56)
        print(N,val,err,'GH difference',abs(val-gh[2]),'reduction bound',red)
