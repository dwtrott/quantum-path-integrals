# Nova, 2026-11-23. Run: python semiclassical_iv.py
import mpmath as mp
from mpmath import iv
iv.dps=40
q=lambda a,b=1:iv.mpf(a)/b
def enclosure(N,T,U,h):
 d=T/N;nu=iv.sqrt(1+U);x=d*nu/2;theta=2*iv.atan2(x,iv.sqrt(1-x*x));assert N*theta<iv.pi
 c0=d/(iv.sin(theta)*iv.sin(N*theta));S2=S3=S4=Q2=Q4=P2=P4=q(0)
 for j in range(1,N):
  l=iv.sin(j*theta);r=iv.sin((N-j)*theta);G=c0*l*r
  S2+=G**2;S3+=G**3;S4+=G**4
  Q2+=G**4+2*c0**2*P2*G*r**2;Q4+=G**4+2*c0**4*P4*r**4
  P2+=G*l**2;P4+=l**4
 S2*=d;S3*=d;S4*=d;Q2*=d*d;Q4*=d*d
 B1=3*U*h*S2/8;B2=h*h*(15*U*S3/48+105*U**2*S2**2/128)
 EA2=U**2*h**2/64*(9*S2**2+72*Q2+24*Q4)
 B3=h**3*(105*U*S4/384+945*U**2*S2*S3/384+10395*U**3*S2**3/3072)
 re=1-15*U*h*h*S3/48-EA2/2;im=-B1
 print('N=',N,'U=',U,'h=',h,'Z2 real=',re,'imag=',im,'B1=',B1,'B2=',B2,'B3=',B3,flush=True)
 return re,im,B1,B2,B3
if __name__=='__main__':
 print(mp.__version__)
 for N,T,U,h in [(2560,q(2),q(3,10),q(1,100)),(2560,q(6,5),q(3),q(1,100)),(2560,q(2),q(1),q(1,1000))]:
  re,im,B1,B2,B3=enclosure(N,T,U,h)
  assert B3<q(2,1000000)
