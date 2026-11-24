# Modification of Nova 028 validated_soft.py, same interval model and weights.
# Run after copying validated_soft.py from its original branch into this directory.
from validated_soft import *

def run_short(N,hs,es,panels,Q='0.3',n=12):
    h=I(hs); eta=I(es); t=time.perf_counter(); _,f,disk,T,a=model(N,h,eta)
    scale=iv.sqrt(iv.sqrt(h)); Q=I(Q); end=Q/scale
    W=weights(n); ww=[rat(x) for x in W]; l=end/(2*panels); R=8*l
    val=iv.mpc(0); err=I(0)
    for j in range(panels):
        c=(2*j+1)*l
        val+=2*l*sum(w*f(c+l*(I(2*k)/n-1)) for k,w in enumerate(ww))
        err+=2*l*(1+sum(rat(abs(w)) for w in W))*disk(c,R)*I(1)/8**(n+2)*I(8)/7
    val*=2; err*=2
    d=T/N; Gnorm=1/((2/d)*(iv.cos((iv.pi+eta)/N)-iv.cos(2*iv.pi/N)))
    M=3*T/4
    c=3*T/128; assert abs(a).b <= (4*c*Q**2).a
    tail=2*(I('0.5')-Q)/scale*iv.exp(abs(a)*Q**2/(2*h)-c*Q**4/h)*(1+(3*T*M+T*d*Gnorm/h)/2)
    print('short',N,es,panels,'value',val,'quadrature',err,'tail',tail,'sec',time.perf_counter()-t,flush=True)
    return val,err,tail
if __name__=='__main__':
    for panels in [12,18,24]: run_short(80,'1e-5','0.001',panels)
