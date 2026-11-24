"""Finite-mesh corrected-overlap enclosure. See accompanying proof for budgets."""
import math
import numpy as np
from classical_scalar import classical

def gaussian(t):
    # Approximate exp(-t), t>=0. Cutoff error below 3e-20.
    z=-np.minimum(t,45.)/1024
    a=np.ones_like(z)
    for k in range(12,0,-1):
        a=1+z*a/k
    for _ in range(10): a=a*a
    return a

def cis(t):
    j=np.rint(t/(np.pi/2)).astype(np.int64)
    s=t-j*(np.pi/2)
    a=np.ones_like(s,dtype=complex)
    for k in range(21,0,-1): a=1+(1j*s)*a/k
    return a*np.array([1,1j,-1,-1j])[j%4]

class Transform:
    def __init__(self,M):
        self.M=M
        ix=np.arange(M); rev=np.zeros(M,dtype=int)
        for _ in range(M.bit_length()-1):
            rev=2*rev+(ix%2);ix//=2
        self.rev=rev
        self.tables=[cis(-np.pi*np.arange(k)/k) for k in (1<<np.arange(M.bit_length()-1))]
    def fft(self,x):
        v=x[self.rev].copy()
        for roots in self.tables:
            k=len(roots); a=v.reshape(-1,2*k)
            y=a[:,k:]*roots; z=a[:,:k].copy()
            a[:,:k]=z+y;a[:,k:]=z-y
        return v
    def ifft(self,x):return self.fft(x.conj()).conj()/self.M

def wave(n,H,transform=None):
    assert H in (20,40,80) and 1<=n<=4096
    M=H//20*4096
    T=Transform(M) if transform is None else transform
    assert T.M==M
    x=(np.arange(M)-M//2)*(48/M)
    k=np.r_[np.arange(M//2),np.arange(-M//2,0)]*(2*np.pi/48)
    potential=3*gaussian(.5*x*x)
    half=cis(-H/n*potential)
    kick=cis(-2*H/n*potential)
    drift=cis(-k*k/(n*H))
    psi=half*gaussian(.5*H*(x+2)**2)*cis(2*H*x)
    for j in range(n):
        psi=T.ifft(drift*T.fft(psi))
        psi*=half if j+1==n else kick
    return x,psi

def overlap(n,H,state,parameters=None,wave_error=0):
    x,psi=state
    M=len(x)
    q,p,act,base,a1=classical(n) if parameters is None else parameters
    mid=lambda v:(v.l+v.h)/2
    assert q.mag()<15 and p.mag()<6.01 and act.mag()*H<2000
    target=gaussian(.5*H*(x-mid(q))**2)*cis(-H*mid(p)*x)
    A=np.sum(target*psi)*(48/M)*np.sqrt(H/np.pi)
    baseval=base.mid()
    L=cis(np.array(H*mid(act)))*baseval
    ratio=A/L-1-a1.mid()/H
    # Normalized overlap and scalar errors: use the actual interval widths.
    errA=wave_error+2e-6+(q.h-q.l)*np.sqrt(H)+(p.h-p.l)*H*25
    errL=1e-10+base.radius()+(act.h-act.l)*H*(abs(baseval)+1)
    assert abs(L)>errL+.1
    error=(errA+abs(A/L)*errL)/(abs(L)-errL)+a1.radius()/H+1e-10
    return abs(ratio)*H*H, (abs(ratio)+error)*H*H
def run(n,H,transform=None):
    return overlap(n,H,wave(n,H,transform))

if __name__=="__main__":
    import sys
    for H in (20,40,80):
        T=Transform(H//20*4096)
        for n in (map(int,sys.argv[1:]) if sys.argv[1:] else (1,2,4,8,16,32,64,128,256,1024)):
            print(n,H,*run(n,H,T),flush=True)
