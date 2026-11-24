"""Test of two candidate selection rules on a grid of targets (floating point, hbar = 0.1 and 0.05).
For each target: all saddles in a box (Newton multistart), the exact amplitude, and the single saddle whose leading term
matches it best.  Rule R0: the dominant saddle is the one with the smallest nonnegative Im Phi.
Rule NR ('near-real'): the dominant saddle has small imaginary parts."""
import numpy as np, sys
from stokes_map import grad, hess, Phi, exact, lead
rng = np.random.default_rng(11)
starts = [rng.uniform(-4, 4, 2) + 1j*rng.uniform(-2.2, 2.2, 2) for _ in range(1500)]
def all_saddles(b, pb):
    found = []
    with np.errstate(all='ignore'):
        for z0 in starts:
            z = z0.copy(); ok = True
            for _ in range(50):
                try: z = z - np.linalg.solve(hess(z), grad(z, b, pb))
                except np.linalg.LinAlgError: ok = False; break
                if not np.all(np.isfinite(z)) or np.max(abs(z)) > 30: ok = False; break
            if ok and np.linalg.norm(grad(z, b, pb)) < 1e-10 and np.max(abs(z.real)) < 6 and np.max(abs(z.imag)) < 4 and all(np.linalg.norm(z-w) > 1e-6 for w in found): found.append(z)
    return found
bs = np.linspace(-7, 5, 13); ps = np.linspace(-5, 4, 10)
stats = dict(n=0, single_ok=0, r0_ok=0, r0_bad=[], imz=[], two=0)
print("target (b,pb) | #saddles Im>=0 | rate from two hbars | dominant Im Phi | smallest nonneg Im Phi | max|Im z| of dominant | fit at hbar=0.05")
for pb in ps:
    for b in bs:
        I1, I2 = exact(b, pb, 0.1), exact(b, pb, 0.05, n=2001)
        if abs(I2) < 1e-11 or abs(I1) < 1e-11: continue
        S = [z for z in all_saddles(b, pb) if Phi(z[0], z[1], b, pb).imag > -1e-6]
        if not S: continue
        S.sort(key=lambda z: Phi(z[0], z[1], b, pb).imag)
        fits = [min(abs(I2/lead(z, b, pb, 0.05) - 1), abs(I2/lead(z, b, pb, 0.05) + 1)) for z in S]
        k = int(np.argmin(fits)); zd = S[k]; imd = Phi(zd[0], zd[1], b, pb).imag; ims = Phi(S[0][0], S[0][1], b, pb).imag
        rate = (np.log(abs(I1)/0.1) - np.log(abs(I2)/0.05))/(1/0.05 - 1/0.1)
        stats['n'] += 1; ok = fits[k] < 0.15
        if ok:
            stats['single_ok'] += 1; stats['imz'].append(np.max(abs(zd.imag)))
            if k == 0: stats['r0_ok'] += 1
            else: stats['r0_bad'].append((b, pb, imd, ims, np.max(abs(S[0].imag)), np.max(abs(zd.imag))))
        else: stats['two'] += 1
        if (stats['n'] % 9) == 1 or (ok and k != 0):
            print(f"({b:+.1f},{pb:+.1f}) | {len(S):3d} | {rate:.4f} | {imd:.4f} | {ims:.4f} | {np.max(abs(zd.imag)):.3f} | {fits[k]:.3f}{'   <-- R0 fails: a saddle with smaller nonnegative Im Phi does not contribute' if ok and k != 0 else ''}")
print(f"\ntargets with resolvable amplitude: {stats['n']};  one saddle reproduces the amplitude within 15% at hbar=0.05: {stats['single_ok']};  not single-saddle (interference or near a turning point): {stats['two']}")
print(f"rule R0 (smallest nonnegative Im Phi) holds in {stats['r0_ok']} of {stats['single_ok']} single-saddle targets; fails in {len(stats['r0_bad'])}")
im = np.array(stats['imz']); print(f"max |Im z| of the dominant saddle: median {np.median(im):.3f}, 90th percentile {np.percentile(im, 90):.3f}, maximum {im.max():.3f}")
for r in stats['r0_bad'][:12]: print(f"   R0 failure at (b,pb)=({r[0]:+.1f},{r[1]:+.1f}): dominant Im Phi {r[2]:.3f}, non-contributing saddle with Im Phi {r[3]:.3f} and max|Im z| {r[4]:.3f} (dominant has max|Im z| {r[5]:.3f})")
