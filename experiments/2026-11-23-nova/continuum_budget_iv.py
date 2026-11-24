"""Run python continuum_budget_iv.py; interval evaluation of continuum error constants, not finite integral."""
from mpmath import iv
iv.dps=40
T=iv.mpf(2); h=1/iv.mpf(10); t=iv.mpf(1); L1=iv.sqrt(2/iv.pi); L2=iv.mpf(1); L3=2*L1; L4=iv.mpf(3)
b1,b2,b3,b4=iv.exp(t)*L1,iv.exp(t)*(L2+L1**2),iv.exp(t)*(L3+3*L1*L2+L1**3),iv.exp(t)*(L4+4*L1*L3+3*L2**2+6*L1**2*L2+L1**4)
for N in [640,2560]:
    d=T/N;theta=2*iv.atan2(d/2,iv.sqrt(1-d*d/4));nu=theta/d;M=iv.sin(theta)/theta;s=iv.sin(T)-T*(nu-1);se=iv.sin(nu*T);G1=1/(M*se);G2=nu/(M*se)
    err=d*d/8*((2*h*G2+h/(M*T))*b2+h*h*G1*G1*b4)+h*b2/2*((1/(M*nu)-1)/s+T*(nu-1)*(1/s+1/s**2))
    assert err < (iv.mpf('2.798e-6') if N==640 else iv.mpf('1.749e-7'))
    print(N,err)
