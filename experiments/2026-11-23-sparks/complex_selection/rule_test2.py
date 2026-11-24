"""Strip rule test (floating point).  For each target, collect all saddles with Im Phi >= 0, find the one that reproduces
the exact amplitude at hbar = 0.05, and test: 'the contributing saddle is the one of smallest nonnegative Im Phi among
saddles with max|Im z_j| < B'.  Saves the raw data to rule_data.pkl."""
import numpy as np, pickle
from rule_test import all_saddles
from stokes_map import Phi, exact, lead
bs = np.linspace(-7, 5, 13); ps = np.linspace(-5, 4, 10); data = []
for pb in ps:
    for b in bs:
        I1, I2 = exact(b, pb, 0.1), exact(b, pb, 0.05, n=2001)
        if abs(I2) < 1e-11 or abs(I1) < 1e-11: continue
        S = [z for z in all_saddles(b, pb) if Phi(z[0], z[1], b, pb).imag > -1e-6]
        with np.errstate(all='ignore'):
            fits = [min(abs(I2/lead(z, b, pb, 0.05) - 1), abs(I2/lead(z, b, pb, 0.05) + 1)) for z in S]
        data.append(dict(b=b, pb=pb, S=S, im=[Phi(z[0], z[1], b, pb).imag for z in S], fits=fits))
pickle.dump(data, open('rule_data.pkl', 'wb'))
single = [d for d in data if d['S'] and min(d['fits']) < 0.15]
print(f"targets resolved {len(data)}, single-saddle {len(single)}")
for B in (0.85, 1.0, 1.2, 1.5, 2.0, 9.9):
    ok = 0; bad = []
    for d in single:
        k = int(np.argmin(d['fits'])); c = [(im, j) for j, (z, im) in enumerate(zip(d['S'], d['im'])) if np.max(abs(z.imag)) < B]
        if c and min(c)[1] == k: ok += 1
        else: bad.append((d['b'], d['pb']))
    print(f"strip half-width B = {B}: rule holds at {ok} of {len(single)}; failures at {bad[:8]}")
# gap between dominant saddle and the nearest excluded saddle, in max|Im z|
gin = max(np.max(abs(d['S'][int(np.argmin(d['fits']))].imag)) for d in single)
gout = min(np.max(abs(z.imag)) for d in single for j, (z, im) in enumerate(zip(d['S'], d['im'])) if im < d['im'][int(np.argmin(d['fits']))] - 1e-9)
print(f"largest max|Im z| over contributing saddles {gin:.3f}; smallest max|Im z| over lower non-contributing saddles {gout:.3f}")
