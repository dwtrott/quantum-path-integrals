"""Certified block penalty constants for the one-mode contour; not global selection."""
import sys,time
from fractions import Fraction as F
import numpy as np
import scipy,mpmath
from mpmath import iv
from structured_audit import setup
from structured_middle_obstruction import root_ball,I
iv.dps=40

def Q(x):
    x=F(x);return I(x.numerator)/x.denominator

def G(d,k,i,j):return d*min(i,j)*(k-max(i,j))/k

def check_exact():
    # Sharp linear-site example and Green inverse / constant resource identity.
    for k in range(2,33):
        d=F(1,7);H=k*d
        g=[[G(d,k,i,j) for j in range(1,k)] for i in range(1,k)]
        for i in range(k-1):
            for j in range(k-1):
                assert (2*g[i][j]-(g[i-1][j] if i else 0)-(g[i+1][j] if i+1<k-1 else 0))/d==int(i==j)
        B=d*d*sum(map(sum,g));assert B==H*(H*H-d*d)/12
        w=[F(0)]+[-sum(row)*d/2 for row in g]+[F(0)]
        E=sum((y-x)**2/d for x,y in zip(w,w[1:]));D=d*sum(w)
        assert E+D==-B/4 and D*D==B*E
    print('exact rational Green and sharpness tests: k=2..32 passed',flush=True)

def f(t,x):return min(F('0.8291')*x,t*x+F('.4')) # x >=0

def certify(N):
    start=time.perf_counter();_,_,_,K,D,_,z,_,_=setup(N);d=F(2,N);dd=Q(d);rho=I('1e-8');R=F('4.8')
    mm=[dd/2 if j in (0,N) else dd for j in range(N+1)]
    zz=[v+iv.mpc(iv.mpf([-1,1])*rho,iv.mpf([-1,1])*rho) for v in root_ball(N,z,K,D,mm,dd,rho)]
    t=[F('.02')+F('1.4')*F(j,N)*(1-F(j,N)) for j in range(N+1)]
    s=F('.8291');B=F('.4');ks=[B/(s-v) for v in t]
    # Potential derivative for all x in [-R,R] and c in [-.46,.46]. Split at the exact kink.
    lp=[I(0)]*(N+1);lc=[I(0)]*(N+1);ls=[I(0)]*(N+1)
    for j in range(1,N):
        Z=zz[j];base=Z*iv.exp(-Z*Z/2);L=I(0)
        for l,r,beta,slope in [(-R,-ks[j],-B,t[j]),(-ks[j],ks[j],F(0),s),(ks[j],R,B,t[j])]:
            for i in range(48):
                a0=l+(r-l)*i/48;a1=l+(r-l)*(i+1)/48
                X=iv.mpf([Q(a0).a,Q(a1).b])
                for k in range(4):
                    C=iv.mpf([Q(F(23,100)*(k-2)).a,Q(F(23,100)*(k-1)).b])
                    W=Z+iv.mpc(X,Q(slope)*X+Q(beta)+C)
                    deriv=(3*dd*(W*iv.exp(-W*W/2)-base)*(1+1j*Q(slope))).imag
                    L=max(L,abs(deriv).b)
        lp[j]=L
        # Exact telescoping site correction derivative is odd, piecewise linear.
        pts=[F(0),R]+[x for x in ks[j-1:j+2] if x<=R]
        lc[j]=Q(max(abs(2*f(t[j],x)-f(t[j-1],x)-f(t[j+1],x))/d for x in pts)).b
        # Optional 042 lower envelope: C_j >= (2t_j-t_{j-1}-t_{j+1})/(2d)*(x²-k_j²)_+.
        ls[j]=Q((2*t[j]-t[j-1]-t[j+1])*R/d).b
    assert I('1.6247')*I('.4')/iv.sqrt(2)<I('.46')
    for stride in [2,4,8,16]:
        if N%stride:continue
        A=I(0);As=I(0);Bm=I(0)
        for a in range(0,N,stride):
            k=stride;t0=Q(min(t[a:a+k+1]));H=k*dd
            for i in range(1,k):
                for j in range(1,k):
                    gij=G(dd,k,i,j)/t0
                    A+=(lp[a+i]+lc[a+i])*(lp[a+j]+lc[a+j])*gij
                    As+=(lp[a+i]+ls[a+i])*(lp[a+j]+ls[a+j])*gij
            Bm+=H*(H*H-dd*dd)/(12*t0)
        assert A>0 and As>0 and Bm>0
        print('N',N,'stride',stride,'nodes',N//stride+1,'exact site loss upper',float((A/4).b),'lower-envelope loss upper',float((As/4).b),'mean coefficient B upper',float(Bm.b),'seconds',time.perf_counter()-start,flush=True)

    # Floating proposal for nonuniform block placement, then interval re-evaluation.
    # The proposal need not be optimal for the certified bound to hold.
    def choose(L,blocks):
        cost=np.full((N+1,N+1),np.inf)
        for a in range(N):
            for b in range(a+1,min(N,a+16)+1):
                k=b-a;v=np.array([float(x.b) for x in L[a+1:b]])
                ii=np.arange(1,k);g=float(d)*np.minimum.outer(ii,ii)*(k-np.maximum.outer(ii,ii))/k
                cost[a,b]=(v@g@v)/(4*float(min(t[a:b+1])))
        dp=np.full((blocks+1,N+1),np.inf);pred=np.full((blocks+1,N+1),-1,dtype=int);dp[0,0]=0
        for j in range(1,blocks+1):
            for b in range(j,N+1):
                a=max(j-1,b-16);vals=dp[j-1,a:b]+cost[a:b,b];i=int(np.argmin(vals));dp[j,b]=vals[i];pred[j,b]=a+i
        ends=[N]
        for j in range(blocks,0,-1):ends.append(int(pred[j,ends[-1]]))
        assert ends[-1]==0 and all(x>y for x,y in zip(ends,ends[1:]))
        ends=ends[::-1];A=I(0);Bm=I(0)
        for a,b in zip(ends,ends[1:]):
            k=b-a;t0=Q(min(t[a:b+1]));H=k*dd
            for i in range(1,k):
                for j in range(1,k):A+=L[a+i]*L[a+j]*G(dd,k,i,j)/t0
            Bm+=H*(H*H-dd*dd)/(12*t0)
        return ends,A/4,Bm
    for stride in [2,4]:
        if N%stride:continue
        for mode,L in [('exact',[x+y for x,y in zip(lp,lc)]),('envelope',[x+y for x,y in zip(lp,ls)])]:
            ends,A,Bm=choose(L,N//stride)
            print('adaptive',N,stride,mode,'loss upper',float(A.b),'mean B upper',float(Bm.b),'ends',ends,flush=True)

if __name__=='__main__':
    print('versions',np.__version__,scipy.__version__,mpmath.__version__,flush=True)
    check_exact()
    for N in map(int,sys.argv[1:] or [16,32,64,128]):certify(N)
