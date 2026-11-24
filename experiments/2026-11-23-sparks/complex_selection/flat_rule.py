"""Flat-plane rule test (floating point).  A saddle z_s = x_s + i y_s is 'flat-selected' when x -> Im Phi(x + i y_s)
attains its global minimum over R^2 at x_s.  Then shifting the integration plane to R^2 + i y_s and Laplace's method give
I ~ leading term of z_s, for every small hbar.  Uses the saddle data of rule_test2.py."""
import numpy as np, pickle
from stokes_map import Phi, hess
data = pickle.load(open('rule_data.pkl', 'rb'))
xs = np.linspace(-9, 9, 721); X1, X2 = np.meshgrid(xs, xs, indexing='ij')
def flat_gap(z, b, pb):
    m = Phi(X1 + 1j*z[0].imag, X2 + 1j*z[1].imag, b, pb).imag.min()
    return Phi(z[0], z[1], b, pb).imag - m          # 0 (up to grid error) when the saddle is the global minimum of its plane
n = dom_flat = 0; counts = []; rows = []
for d in data:
    if not d['S'] or min(d['fits']) >= 0.15: continue
    n += 1; k = int(np.argmin(d['fits'])); b, pb = d['b'], d['pb']
    gaps = [flat_gap(z, b, pb) for z in d['S']]; flat = [j for j, gp in enumerate(gaps) if gp < 1e-3]
    counts.append(len(flat)); ok = (flat == [k]); dom_flat += ok
    ev = np.linalg.eigvalsh(hess(d['S'][k]).imag)
    rows.append((b, pb, d['im'][k], gaps[k], ev.min(), len(flat), ok))
print(f"single-saddle targets: {n};  the contributing saddle is the unique flat-selected saddle at {dom_flat} of them")
print(f"number of flat-selected saddles per target: min {min(counts)}, max {max(counts)}")
print("targets where the flat rule does not decide (b, pb, Im Phi of contributing saddle, depth of plane minimum below it, min eig of Im Hess):")
for r in rows:
    if not r[6]: print(f"   ({r[0]:+.1f},{r[1]:+.1f})  Im Phi_s = {r[2]:.4f}  gap = {r[3]:.4f}  min eig Im Phi'' = {r[4]:+.3f}  flat-selected saddles: {r[5]}")
