"""Reproduce one tranche of the direct spectral certificate, with assertions."""
import sys
from spectral_corrected import run, Transform

def check(lo,hi):
    assert 1 <= lo <= hi <= 512
    for H in (20,40,80):
        transform=Transform(H//20*4096)
        peak=max((run(n,H,transform)[1],n) for n in range(lo,hi+1))
        print(H,peak,flush=True)
        assert peak[0]<1.6

if __name__=="__main__":
    check(*map(int,sys.argv[1:])) if sys.argv[1:] else check(1,512)
