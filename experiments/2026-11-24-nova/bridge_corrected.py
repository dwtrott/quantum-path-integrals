from spectral_corrected import wave,overlap
from classical_scalar import classical
from uniform_049 import I
def bridge(lo=513,hi=1023):
    ref=4096
    data={h:wave(ref,h) for h in (20,40,80)}
    peaks={h:(0,0) for h in data}
    for n in range(lo,hi+1):
        params=classical(n)
        for H,state in data.items():
            err=4*H*((I(2)/n).sq()+(I(2)/ref).sq())
            ub=overlap(n,H,state,params,err.h)[1]
            peaks[H]=max(peaks[H],(ub,n))
    print(peaks,flush=True)
if __name__=='__main__': bridge()
