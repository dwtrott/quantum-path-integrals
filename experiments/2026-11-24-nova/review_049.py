import sys
import numpy as np

N, r, g, T = 64, .2, 3., 2.
d = T/N
m = np.full(N+1, d)
m[[0,-1]] /= 2
x = np.empty(N+1)
x[0], p = -2., 2. + g*m[0]*(-2.)*np.exp(-2.)
for j in range(1,N+1):
    x[j] = x[j-1] + d*p
    p += g*m[j]*x[j]*np.exp(-x[j]**2/2)
def derivatives(z):
    W = np.exp(-z*z/2)
    return (3*z-z**3)*W, (z**4-6*z*z+3)*W
D = np.diff(np.eye(N+1),axis=0)
H = (D.T@D/d - np.diag(g*m*(x*x-1)*np.exp(-x*x/2))).astype(complex)
H[0,0] += 1j
H[-1,-1] += 1j
A = -1j*H
B, C = A.real, A.imag
lam, V = np.linalg.eigh(np.block([[B,C],[C,-B]]))
sigma = lam[-(N+1):]
U = V[:N+1,-(N+1):] - 1j*V[N+1:,-(N+1):]
R = U/np.sqrt(sigma)
P = (R.conj().T @ (m[:,None]*R)).real
mu, V = np.linalg.eigh(P)
s2 = np.sum(abs(R)**2,axis=1)
G = np.linalg.solve(H,np.eye(N+1))
gd = np.diag(G)
w3,w4 = np.array(derivatives(x))*m
a1 = 1j*g/8*(w4@gd**2) + 1j*g*g/72*(9*(w3*gd)@G@(w3*gd)+6*w3@(G**3)@w3)
res = D.T@(D@x)/d + g*m*x*np.exp(-x*x/2)
res[0] += 2
res[-1] -= p
print(sys.version, np.__version__)
print('qb pb stationarity Takagi orthogonality',x[-1],p,max(abs(res)),np.linalg.norm(R.T@A@R-np.eye(N+1),2),np.linalg.norm(U.conj().T@U-np.eye(N+1),2))
for nx,ny in [(2401,41),(12001,81)]:
    z = np.linspace(-6,6,nx)[None,:] + 1j*np.linspace(-r,r,ny)[:,None]
    M3,M4 = [np.max(abs(v)) for v in derivatives(z)]
    k = g*M3*r/3
    F = np.exp(-.5*np.log1p(-k*mu).sum())
    sp2 = (abs(R@V)**2)@(1/(1-k*mu))
    C1 = g*M4/8*(m@s2**2)
    C2 = F*(g*M3)**2*15/72*(m@sp2**1.5)**2
    print(nx,ny,'M3 M4 normP traceP maxs2 kNorm F C1 C2 sum |a1|')
    print(M3,M4,mu[-1],mu.sum(),s2.max(),k*mu[-1],F,C1,C2,C1+C2,abs(a1))


# Independently compute the principal matrix inverse root by eigendecomposition.
e, Q = np.linalg.eig(A)
R = (Q*(1/np.sqrt(e)))@np.linalg.solve(Q,np.eye(N+1))
P = (R.conj().T @ (m[:,None]*R)).real
mu, V = np.linalg.eigh(P)
s2 = np.sum(abs(R)**2,axis=1)
print('principal-root eigenvector condition / congruence residual',
      np.linalg.cond(Q), np.linalg.norm(R.T@A@R-np.eye(N+1),2))
for nx,ny in [(2401,41),(12001,81)]:
    z = np.linspace(-6,6,nx)[None,:] + 1j*np.linspace(-r,r,ny)[:,None]
    M3,M4 = [np.max(abs(v)) for v in derivatives(z)]
    k = g*M3*r/3
    F = np.exp(-.5*np.log1p(-k*mu).sum())
    sp2 = (abs(R@V)**2)@(1/(1-k*mu))
    C1 = g*M4/8*(m@s2**2)
    C2 = F*(g*M3)**2*15/72*(m@sp2**1.5)**2
    print(nx,ny,'M3 M4 normP traceP maxs2 kNorm F C1 C2 sum |a1|')
    print(M3,M4,mu[-1],mu.sum(),s2.max(),k*mu[-1],F,C1,C2,C1+C2,abs(a1))
