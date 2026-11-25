"""Independent floating point evaluation of the constant C_N of Nova 050 (unitary Duhamel bound at the real classical target).
C_N = sqrt(|w|/2) (9 S4/32 + 345 S6/576),  S4 = sum m_j |L_j|^2 |R_j|^2,  S6 = sum m_j^2 |L_j|^3 |R_j|^3 + 2 sum_{i<j} m_i m_j |L_i|^3 |R_j|^3,
L, R: position components of the forward and reversed linearised Verlet flow started from (1, i).  Also |A0| = sqrt(2/|w|) against the
oriented leading term of my own code, and the bound C_N hbar against the measured |I/L - 1|."""
import numpy as np
from nslice import g, setup, classical, W, Ttot
from order_hbar import a1, saddle
V2 = lambda x: g*(x*x-1)*W(x)
for N in (16, 32, 64, 128, 256, 1024):
    d, w, K = setup(N); m = d*w; x, p = classical(N); c = -V2(x)
    y, v = 1+0j, 1j; L = np.zeros(N+1, complex)
    for j in range(N+1):
        L[j] = y; v = v + m[j]*c[j]*y
        if j < N: y = y + d*v
    wr = v + 1j*y
    y, v = 1+0j, 1j; R = np.zeros(N+1, complex)
    for j in range(N, -1, -1):
        R[j] = y; v = v + m[j]*c[j]*y
        if j > 0: y = y + d*v
    aL, aR = abs(L), abs(R); S4 = np.sum(m*aL**2*aR**2); cum = np.cumsum(m*aL**3)
    S6 = np.sum(m*m*aL**3*aR**3) + 2*np.sum((cum[:-1])*(m[1:]*aR[1:]**3))
    C = np.sqrt(abs(wr)/2)*(9*S4/32 + 345*S6/576)
    extra = ""
    if N <= 256:
        z, qb, pb = saddle(N, 0, 0); extra = f"  |a_1| = {abs(a1(z, N)):.4f}"
    print(f"N = {N:5d}: |w| = {abs(wr):.5f}, S4 = {S4:.4f}, S6 = {S6:.4f}, C_N = {C:.4f}{extra}")
