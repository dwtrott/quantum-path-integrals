from nova_transfer import *

def pa_lattice(N,s=.25,M=32,L=14):
    d,lam,_=params(N); n=N-1; tau=U*d/h; t=n*abs(tau); r=2
    j=np.arange(1,N); v=np.sqrt(2/N)*np.sin(np.pi*j[:,None]*np.arange(1,r+1)[None,:]/N)
    a=np.zeros(N); a[r+1:]=1/lam[r:]
    diag=(np.sum(a)-np.fft.fft(a).real[1:N])/N
    D=1+1j*h*diag; b=np.sqrt(2*h)*v/np.sqrt(lam[:r]); zfactor=np.exp(1j*np.pi/4)*b
    # Safe spectral bounds, not rounded upwards.
    B=2*h/N*np.sum(1/lam[:r]); K=h*t*t*B*2/N*np.sum(1/lam[r:])
    xs=s*np.arange(M+1); ws=np.full(M+1,2.); ws[0]=1.
    inv=np.zeros(n); inv[r:]=1/lam[r:]
    result=0j
    for x,wx in zip(xs,ws):
        q=np.column_stack((np.full(len(xs),x),xs)); z=q@zfactor.T
        W=np.exp(-z*z/(2*D))/np.sqrt(D); wp=-z*W/D
        F=np.exp(1j*tau*W.sum(axis=1)); wd=dst(wp,type=1,norm='ortho',axis=1)
        P=1j*h*(1j*tau)**2*np.sum(wd*wd*inv,axis=1)/2
        result+=wx*np.dot(ws,F*(1+P)*np.exp(-np.sum(q*q,axis=1)))
    a0=np.pi**2/(s*s*(1+(L+2)*B)**2); v0=np.exp(-a0); ltheta=np.log1p(2*v0/(1-v0)); theta=np.exp(ltheta)
    alias=np.exp(t)*((1+K/2)*np.expm1(r*ltheta)+K*np.pi**2/s**2*r*theta**(r-1)*2*v0*(1+v0)/(1-v0)**3)
    R=t**(L+1)/math.factorial(L+1)/(1-t/(L+2))
    tay=R*(1+K*r/2+(1+2*K/np.e)*(s+np.sqrt(2*np.pi))**r/np.pi**(r/2))
    box=np.exp(t)*(1+2*K/np.e)*r*2*np.exp(-s*s*M*M/2)/(s*M)*(s+np.sqrt(2*np.pi))**(r-1)/np.pi**(r/2)
    red=h*h*tau*tau*np.exp(t)*np.sum(inv*inv)/4
    return result*s*s/np.pi,red+alias+tay+box

if __name__=='__main__':
    from statistics import median
    for N in (8,40,160,640,2560):
        row=[]; vals=[]
        for f in (lambda:transfer(N,eps=1e-5),lambda:transfer(N,eps=1e-8),lambda:pa_lattice(N)):
            times=[]
            for rep in range(5):
                t0=time.perf_counter(); q=f(); times.append(time.perf_counter()-t0)
            row.append(median(times)); vals.append(q)
        print(N,'medians',row,'matched grids',vals[0][2],vals[1][2],'bounds',*[v[1] for v in vals], 'diff matched-ref',abs(vals[0][0]-vals[1][0]),'diff PA-ref',abs(vals[2][0]-vals[1][0]),flush=True)
