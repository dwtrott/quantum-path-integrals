from poly_soft import *

def check(N,eta):
    h=I('1e-5'); eta=I(eta); f,pe,disk,T,a,g,G,v=model(N,h,eta,6)
    alpha=iv.pi/N; theta=(iv.pi+eta)/N; d=T/N; tau=d/h
    # Independent spectral contractions, no rank-two generator formula.
    S=[[iv.sin(j*k*alpha) for j in range(1,N)] for k in range(2,N)]
    L=[(2/d)*(iv.cos(theta)-iv.cos(k*alpha)) for k in range(2,N)]
    diag=[sum(2*h/N*row[j]**2/lam for row,lam in zip(S,L)) for j in range(N-1)]
    scale=iv.sqrt(iv.sqrt(h)); rot=iv.exp(-J*iv.pi/8)
    for q0 in ['0.05','0.2','0.3']:
        x=I(q0)/scale; q=I(q0)*rot; r=[q*y for y in v]; z=[-J*y*y/2 for y in r]; E=[iv.exp(y) for y in z]
        w=J*tau*sum(e-1-y for e,y in zip(E,z)); wp=[tau*y*(e-1) for y,e in zip(r,E)]; wpp=[tau*(e-1-J*y*y*e) for y,e in zip(r,E)]
        quad=sum(2*h/N*sum(y*u for y,u in zip(row,wp))**2/lam for row,lam in zip(S,L))
        exact=rot*iv.exp(-a*q*q/(2*h)+w)*(1+(sum(x*y for x,y in zip(wpp,diag))+quad)/2)
        error=abs(f(x)-exact); bound=pe(x)
        assert error.b<=bound.a
        print('independent',N,q0,'error',error,'bound',bound,flush=True)
if __name__=='__main__':
    check(4,'0'); check(80,'0.001')
