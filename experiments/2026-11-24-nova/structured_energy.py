"""Interval envelope for the finite-dimensional one-mode contour. OPENBLAS_NUM_THREADS=1."""
import sys,time
import numpy as np
import scipy,mpmath
from mpmath import iv
from structured_audit import setup
from structured_middle_obstruction import root_ball,I
iv.dps=40

def certify(N):
    start=time.perf_counter();d,m,E,K,D,t,z,_,_=setup(N);n=N+1;dd=I(2)/N;rho=I('1e-8')
    mm=[dd/2 if j in (0,N) else dd for j in range(n)]
    tt=[I('.02')+I('1.4')*(I(j)/N)*(1-I(j)/N) for j in range(n)]
    zz=[v+iv.mpc(iv.mpf([-1,1])*rho,iv.mpf([-1,1])*rho) for v in root_ball(N,z,K,D,mm,dd,rho)]
    sig=I('.8291');B=I('.4');alpha=I('.1');cmax=I('1.6247')*B/iv.sqrt(2);beta=B+cmax
    W=lambda w:iv.exp(-w*w/2)
    C=sum((2*tt[j]-tt[j-1]-tt[j+1])*(B/(sig-tt[j]))**2/(2*dd) for j in range(1,N))
    av=[(tt[j]+tt[j+1])/2 for j in range(N)]
    diag=[2*tt[j]*(1 if j in (0,N) else 2)/dd+(1-tt[j]**2 if j in (0,N) else 0) for j in range(n)]
    off=[-2*v/dd for v in av]
    def ldl(shift):
        piv=[diag[0]-shift*mm[0]];assert piv[0]>0
        for j in range(1,n):
            piv.append(diag[j]-shift*mm[j]-off[j-1]**2/piv[-1]);assert piv[-1]>0
        return piv
    ldl(I('.89'));piv=ldl(2*alpha)
    inv=[1/piv[-1]]
    for j in range(N-1,-1,-1):inv.append(1/piv[j]+(off[j]/piv[j])**2*inv[-1])
    inv=inv[::-1]
    # Solve B_energy u=P^T via interval LDL and compute conditional variances.
    P=[v/iv.sqrt(2) for v in mm];f=[P[0]]
    for j in range(1,n):f.append(P[j]-off[j-1]/piv[j-1]*f[-1])
    u=[f[-1]/piv[-1]]
    for j in range(N-1,-1,-1):u.append((f[j]-off[j]*u[-1])/piv[j])
    u=u[::-1];den=sum(P[j]*u[j] for j in range(n));assert den>0
    var=[inv[j]-u[j]**2/den for j in range(n)];assert all(v>0 for v in var)
    count=0;total=C;maxdepth=0
    for j in range(n):
        endpoint=j in (0,N);g=-3*mm[j]*zz[j]*W(zz[j]);ell=abs(g.imag+tt[j]*g.real)+(beta*tt[j] if endpoint else 0)
        e0=beta*abs(g.real)+3*mm[j]*(iv.exp((abs(zz[j].imag-tt[j]*zz[j].real)+beta)**2/(2*(1-tt[j]**2)))+abs(W(zz[j])))+(beta**2/2 if endpoint else 0)
        R=8
        while not (alpha*mm[j]*R*R-ell*R-e0>0 and 2*alpha*mm[j]*R>ell):R*=2
        x=np.linspace(-R,R,2001)[:,None];c=np.linspace(-.46,.46,81)[None,:]
        b=np.clip((.8291-t[j])*x,-.4,.4)+c;y=t[j]*x+b;w=x+1j*y
        v=3*m[j]*(np.exp(-z[j]**2/2)-np.exp(-(z[j]+w)**2/2)-z[j]*np.exp(-z[j]**2/2)*w).imag+.1*m[j]*x*x
        if endpoint:v-=t[j]*x*b+b*b/2
        bound=int(np.floor((min(0.,float(v.min()))-.4*m[j]-(.04 if endpoint else 0))*10**6));target=I(bound)/10**6
        wz=W(zz[j]);stack=[(-R,R,-.46,.46,0)];assert cmax<iv.mpf(.46)
        while stack:
            l,r,a,b,dep=stack.pop();count+=1;maxdepth=max(maxdepth,dep)
            xx=iv.mpf([l,r]);cc=iv.mpf([a,b]);raw=(sig-tt[j])*xx
            # Enclose the monotone clip using directed endpoints, including uncertain B.
            bb=iv.mpf([max(-B.b,min(B.a,raw.a)),max(-B.a,min(B.b,raw.b))])+cc
            yy=tt[j]*xx+bb;ww=iv.mpc(xx,yy)
            val=3*mm[j]*(wz-W(zz[j]+ww)-zz[j]*wz*ww).imag+alpha*mm[j]*xx**2
            if endpoint:val-=tt[j]*xx*bb+bb**2/2
            if val>=target:continue
            if dep>=45:raise RuntimeError((N,j,l,r,a,b,dep,val,target))
            if r-l>=b-a:
                mid=(l+r)/2;stack.extend([(l,mid,a,b,dep+1),(mid,r,a,b,dep+1)])
            else:
                mid=(a+b)/2;stack.extend([(l,r,a,mid,dep+1),(l,r,mid,b,dep+1)])
        total-=target
    assert total<I('7.24'),total
    radius=iv.sqrt(2*total/(I('.89')-2*alpha));coord=iv.sqrt(2*total*max(x.b for x in inv))
    qzero=iv.sqrt(2*total*max(v.b for v in var));qbound=iv.sqrt(2*total*den)
    assert radius<I('4.6') and coord<I('4.8'),(radius,coord)
    print('conditional q=0 coordinate upper',float(qzero.b),'q upper',float(qbound.b),flush=True)
    print('N',N,'panels',count,'maxdepth',maxdepth,'C',float(C.b),'total upper',float(total.b),'radius upper',float(radius.b),'coordinate upper',float(coord.b),'seconds',time.perf_counter()-start,flush=True)

if __name__=='__main__':
    print('versions',np.__version__,scipy.__version__,mpmath.__version__,flush=True)
    for N in map(int,sys.argv[1:] or [16,32,64,128]):certify(N)
