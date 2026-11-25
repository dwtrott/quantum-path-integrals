"""Exact rational check of ordered residual identity, arbitrary inexact states."""
from fractions import Fraction as F
from random import Random
rng=Random(52)
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def mv(A,a):return tuple(dot(row,a) for row in A)
def tr(A):return tuple(zip(*A))
def vec():return tuple(F(rng.randrange(-8,9),7) for _ in range(2))
def rot():
    t=F(rng.randrange(-4,5),5);c=(1-t*t)/(1+t*t);s=2*t/(1+t*t)
    return ((c,-s),(s,c))
for _ in range(100):
    m=rng.randrange(1,9);U=[rot() for j in range(m)]
    u,v=vec(),vec();a=[vec() for j in range(m+1)];b=[vec() for j in range(m+1)]
    r=[sub(u,a[0])]+[sub(mv(U[j],a[j]),a[j+1]) for j in range(m)]
    s=[sub(mv(tr(U[j]),b[j+1]),b[j]) for j in range(m)]+[sub(v,b[m])]
    def propagate(x,start,stop):
        for j in range(start,stop):x=mv(U[j],x)
        return x
    exact=dot(v,propagate(u,0,m))
    local=dot(b[m],a[m])+dot(s[m],a[m])+sum(dot(b[i],r[i]) for i in range(m+1))
    pair=sum(dot(s[j],propagate(r[i],i,j)) for i in range(m+1) for j in range(i,m+1))
    assert exact==local+pair
print('100 exact rational identity checks passed')
