from large_corrected import all_large
from classical_scalar import classical,Z
from uniform_049 import I,sqrt_interval
def nonnegative_root(x):return sqrt_interval(I(max(0,x.l),x.h))
def magnitude(z):return nonnegative_root(z.r.sq()+z.i.sq())
def family(n0=512,slabs=4096,ref=4096):
    q,p,L,a,action=all_large(n0,slabs)
    q0,p0,_,L0,a0=classical(ref)
    dx,dp=q-q0,p-p0
    ar=all_large(ref,slabs)[-1]
    # Both discretizations compared to the same smooth extremal.
    # pbar-p_cont is bounded by the entire terminal p interval width.
    phase=I(action)+ar+(I(p.h)-p.l)*dx.mag()
    for H in (20,40,80):
        coherent=nonnegative_root((dx.sq()+dp.sq())*H/2)+H*phase
        quantum=4*H*((I(2)/n0).sq()+(I(2)/ref).sq())
        diff=L0*(Z(1)+a0/H)-L*(Z(1)+a/H)
        eps=(quantum+coherent+magnitude(diff)+magnitude(L0)*1.6/(H*H))/magnitude(L).l
        print(n0,H,'coherent',coherent.h,'Strang',quantum.h,'constant',(eps*(H*H)).h,'error',eps.h,flush=True)
if __name__=='__main__':
    import sys
    family(*map(int,sys.argv[1:]))
