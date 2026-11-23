import numpy as np, mpmath as mp
mp.mp.dps = 30
# (1) discrete harmonic oscillator: det of interior Hessian and discrete caustics
m=1.0; w=1.0; hbar=1.0
def discrete_K(a,b,T,N):
    dt=T/N
    n=N-1
    H=(m/dt)*(np.diag((2-w*w*dt*dt)*np.ones(n))-np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1))
    # S = 1/2 x^T H x - c^T x + S0 ; c from endpoints
    c=np.zeros(n); c[0]+= (m/dt)*a; c[-1]+=(m/dt)*b
    S0=(m/(2*dt))*(a*a+b*b)-dt*0.5*m*w*w*a*a
    xs=np.linalg.solve(H,c)
    Scl=0.5*xs@H@xs-c@xs+S0
    ev=np.linalg.eigvalsh(H); neg=(ev<0).sum()
    theta=np.arccos(1-w*w*dt*dt/2)
    det_formula=np.sin(N*theta)/np.sin(theta)
    det_num=np.prod(ev*dt/m)
    # K_N = (m/(2 pi i hbar dt))^{N/2} (2 pi i hbar)^{n/2} |det|^{-1/2} e^{-i pi neg/2} e^{i Scl}
    pref=np.sqrt(m/(2*np.pi*hbar*dt))*np.exp(-1j*np.pi/4)*np.sqrt(abs(np.sin(theta)/np.sin(N*theta)))*np.exp(-1j*np.pi*neg/2)
    return pref*np.exp(1j*Scl/hbar), det_formula/det_num, neg
def exact_K(a,b,T):
    s=np.sin(w*T); k=int(np.floor(w*T/np.pi))
    return np.sqrt(m*w/(2*np.pi*hbar*abs(s)))*np.exp(-1j*np.pi/4-1j*np.pi*k/2)*np.exp(1j*m*w*((a*a+b*b)*np.cos(w*T)-2*a*b)/(2*hbar*s))
a,b=0.3,-0.7
for T in (1.0,4.0):
    Ke=exact_K(a,b,T)
    for N in (8,32,128,512):
        K,r,neg=discrete_K(a,b,T,N)
        print(f"T={T} N={N:5d} det ratio={r:.12f} neg={neg} |K_N-K|={abs(K-Ke):.3e}  N*err={N*abs(K-Ke):.4f}")
# (2) one-variable quartic slice: I = int exp(i(alpha x^2 - beta x^4)) dx, three real saddles
def I_num(alpha,beta):
    ph=mp.e**(-1j*mp.pi/8)
    f=lambda y: mp.e**(1j*(alpha*(ph*y)**2-beta*(ph*y)**4))*ph
    return mp.quad(f,mp.linspace(-12/beta**0.25,12/beta**0.25,400))
def I_saddle(alpha,beta):
    c=mp.sqrt(mp.pi/alpha)*mp.e**(1j*mp.pi/4)
    o=2*mp.sqrt(mp.pi/(2*alpha))*mp.e**(-1j*mp.pi/4)*mp.e**(1j*alpha**2/(4*beta))
    return c,o
for beta in (0.05,0.02,0.01,0.005):
    alpha=1.0
    I=I_num(alpha,beta); c,o=I_saddle(alpha,beta)
    print(f"beta={beta}: I={mp.nstr(I,8)}  central={mp.nstr(c,6)} outer={mp.nstr(o,6)}  |I-central|/|I|={mp.nstr(abs(I-c)/abs(I),4)}  |I-(central+outer)|/|I|={mp.nstr(abs(I-c-o)/abs(I),4)}")
