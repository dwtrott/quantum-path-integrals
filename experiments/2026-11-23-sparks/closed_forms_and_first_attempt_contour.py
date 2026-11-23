import numpy as np, mpmath as mp
mp.mp.dps=30
m=w=hbar=1.0
# (A) closed forms for the discrete classical action, left-point and trapezoid rules; spectrum; discrete caustics
def data(a,b,T,N,rule):
    dt=T/N; n=N-1
    A=(m/dt)*(np.diag((2-w*w*dt*dt)*np.ones(n))-np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1))
    c=np.zeros(n); c[0]+=(m/dt)*a; c[-1]+=(m/dt)*b
    S0=(m/(2*dt))*(a*a+b*b)-(dt*0.5*m*w*w*a*a if rule=='left' else dt*0.25*m*w*w*(a*a+b*b))
    xs=np.linalg.solve(A,c); S=0.5*xs@A@xs-c@xs+S0
    th=np.arccos(1-w*w*dt*dt/2)
    Sf=m*np.sin(th)/(2*dt*np.sin(N*th))*((a*a+b*b)*np.cos(N*th)-2*a*b)+(m*w*w*dt/4*(b*b-a*a) if rule=='left' else 0)
    ev=np.sort(np.linalg.eigvalsh(A)); evf=np.sort((2*m/dt)*(np.cos(th)-np.cos(np.arange(1,N)*np.pi/N)))
    xf=(a*np.sin((N-np.arange(1,N))*th)+b*np.sin(np.arange(1,N)*th))/np.sin(N*th)
    G=np.linalg.inv(A); j=np.arange(1,N); Gf=(dt/m)*np.sin(j*th)*np.sin((N-j)*th)/(np.sin(th)*np.sin(N*th))
    return abs(S-Sf), np.max(abs(ev-evf))/np.max(abs(ev)), np.max(abs(xs-xf)), (ev<0).sum(), int(np.floor(N*th/np.pi)), np.max(abs(np.diag(G)-Gf))
for T in (1.0,4.0,7.5):
    for rule in ('left','trap'):
        print(T,rule,["%.2e"%v if isinstance(v,float) else v for v in data(0.3,-0.7,T,40,rule)])
N=40
print("discrete caustic times w*T=2N sin(k pi/2N), k=1,2:",[2*N*np.sin(k*np.pi/(2*N)) for k in (1,2)],"continuum:",np.pi,2*np.pi)
# (B) N=2 check of Fourier-expansion bound and first-order formula, Gaussian W
V0=1.0
def I_full(eps,sig,a,b,T):
    dt=T/2; lam=(m/dt)*(2-w*w*dt*dt); c=(m/dt)*(a+b); xs=c/lam; tau=eps*dt/hbar
    ph=mp.e**(1j*mp.pi/4)
    f=lambda y: mp.e**(1j*(lam*(xs+ph*y)**2/2-c*(xs+ph*y))/hbar-1j*tau*V0*mp.e**(-(xs+ph*y)**2/(2*sig**2)))*ph
    I=mp.quad(f,mp.linspace(-14*mp.sqrt(hbar/lam),14*mp.sqrt(hbar/lam),120))
    I0=mp.sqrt(2j*mp.pi*hbar/lam)*mp.e**(-1j*c*c/(2*lam*hbar))
    s2=sig**2+1j*hbar/lam
    I1=-1j*tau*I0*V0*sig/mp.sqrt(s2)*mp.e**(-xs**2/(2*s2))
    return I,I0,I1,tau
for sig in (0.3,1.0):
  for eps in (0.01,0.5,2.0):
    I,I0,I1,tau=I_full(eps,sig,0.3,-0.7,1.0)
    print(f"sig={sig} eps={eps}: |I-I0|/|I0|={mp.nstr(abs(I-I0)/abs(I0),6)}  bound e^(tau V0)-1={mp.nstr(mp.e**(tau*V0)-1,6)}  |I1|/|I0|={mp.nstr(abs(I1/I0),6)}  |I-I0-I1|/|I0|={mp.nstr(abs(I-I0-I1)/abs(I0),4)}  2nd-order tail bound={mp.nstr(mp.e**(tau*V0)-1-tau*V0,4)}")
