"""Tilted-plane rule test (floating point).  Surface z = z_s + (1 + i t_j) u_j, u real, T = diag(t1, t2), |t_j| < 1
(so every coordinate stays in the sector |Im| < |Re| at infinity, where the kicks exp(-z^2/2) stay bounded).
A saddle is tilt-selected by T when h(u) = Im Phi(z_s + (1+iT)u) has its only global minimum at u = 0.
Tested on the single-saddle targets of rule_test2.py where the flat rule (T = 0) is silent."""
import numpy as np, pickle
from stokes_map import Phi, hess
data = pickle.load(open('rule_data.pkl', 'rb'))
us = np.linspace(-12, 12, 321); U1, U2 = np.meshgrid(us, us, indexing='ij')
uf = np.linspace(-12, 12, 1201); F1, F2 = np.meshgrid(uf, uf, indexing='ij')
def gap(z, b, pb, t1, t2, fine=False):
    A, B = (F1, F2) if fine else (U1, U2)
    with np.errstate(all='ignore'):
        h = Phi(z[0] + (1+1j*t1)*A, z[1] + (1+1j*t2)*B, b, pb).imag
    return Phi(z[0], z[1], b, pb).imag - np.nanmin(h)
ts = np.linspace(-0.9, 0.9, 19)
single = [d for d in data if d['S'] and min(d['fits']) < 0.15]
und = 0; solved = 0; false_sel = 0
print("target | Im Phi_s | flat gap | best tilt (t1,t2) | gap at best tilt (fine grid) | tilt-selected | other saddles tilt-selected by any T on the grid")
for d in single:
    k = int(np.argmin(d['fits'])); b, pb = d['b'], d['pb']; z = d['S'][k]
    g0 = gap(z, b, pb, 0, 0)
    if g0 < 1e-3: continue
    und += 1
    G = np.array([[gap(z, b, pb, t1, t2) for t2 in ts] for t1 in ts]); i, j = np.unravel_index(np.argmin(G), G.shape)
    # local refinement of the tilt
    best = (G[i, j], ts[i], ts[j])
    for t1 in np.linspace(ts[i]-0.1, ts[i]+0.1, 9):
        for t2 in np.linspace(ts[j]-0.1, ts[j]+0.1, 9):
            if abs(t1) < 1 and abs(t2) < 1:
                gg = gap(z, b, pb, t1, t2)
                if gg < best[0]: best = (gg, t1, t2)
    gf = gap(z, b, pb, best[1], best[2], fine=True); ok = gf < 1e-3; solved += ok
    # false selections: any saddle with smaller Im Phi selected by some tilt on the coarse grid?
    fs = 0
    for jj, w in enumerate(d['S']):
        if jj != k and d['im'][jj] < d['im'][k] and np.max(abs(w.imag)) < 2.5:
            if min(gap(w, b, pb, t1, t2) for t1 in ts[::3] for t2 in ts[::3]) < 1e-3: fs += 1
    false_sel += fs
    print(f"({b:+.0f},{pb:+.0f}) | {d['im'][k]:.4f} | {g0:.4f} | ({best[1]:+.3f},{best[2]:+.3f}) | {gf:.5f} | {'yes' if ok else 'no'} | {fs}")
print(f"\nflat rule silent at {und} targets; a diagonal tilt selects the contributing saddle at {solved} of them; false selections {false_sel}")
