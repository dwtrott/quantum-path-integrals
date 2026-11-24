"""Small N interval transfer benchmark; independent of FFT.
Run: python nova_transfer_iv.py
"""
import math,time,mpmath as mp
iv=mp.iv; iv.dps=30
N=8; n=N-1; h=iv.mpf(1)/10; d=iv.mpf(2)/N; tau=iv.mpf(1)/N
s=iv.mpf(1)/16; M=36; L=18
lam=[d*(4/d**2*iv.sin(iv.pi*k/(2*N))**2-1) for k in range(1,N)]
vmax=2*h/N*sum(1/x for x in lam); Lambda=4/(d*h)+L**2*vmax
a=2*iv.pi**2/(s*s*Lambda); v=iv.exp(-a); theta=1+2*v/(1-v); t=n*tau
alias=iv.exp(t)*(theta**n-1)
R=t**(L+1)/math.factorial(L+1)/(1-t/(L+2)); tay=R*(1+theta**n)
box=iv.exp(t)*theta**n*2*n*iv.exp(-(s*M)**2/(2*vmax))
c=1/(2*h*d); pref=iv.sqrt(c/iv.pi)
grid=[s*j for j in range(-M,M+1)]
weights=[iv.exp(d*y*y/(2*h)+iv.j*tau*iv.exp(-iv.j*y*y/2)) for y in grid]
g=[pref*iv.exp(-c*y*y) for y in grid]
kernel=[s*pref*iv.exp(-c*(s*j)**2) for j in range(-2*M,2*M+1)]
f=[x*w for x,w in zip(g,weights)]
for j in range(1,n):
    f=[w*sum(kernel[2*M+i-k]*fk for k,fk in enumerate(f)) for i,w in enumerate(weights)]
den=pref/iv.sqrt(math.prod(x*d for x in lam))
z=s*sum(x*y for x,y in zip(g,f))/den
print(mp.__version__,'N',N,'M',M,'step',s,'bound',alias+tay+box,'lattice',z)
