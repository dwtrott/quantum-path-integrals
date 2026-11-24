"""Cross-check the polynomial moment identity used in the exterior proof."""
import numpy as np

def check(n):
    d = 2 / n
    mass = np.full(n + 1, d)
    mass[[0, -1]] /= 2
    x, momentum = -2., 2.
    transfer = np.eye(2)
    rows, curvature = [], []
    for j, step in enumerate(mass):
        rows.append(transfer[0].copy())
        c = 3 * (1 - x*x) * np.exp(-x*x / 2)
        curvature.append(c)
        transfer = np.array([[1, 0], [step*c, 1]]) @ transfer
        momentum += step * 3*x*np.exp(-x*x / 2)
        if j < n:
            x += d * momentum
            transfer = np.array([[1, d], [0, 1]]) @ transfer

    a, b, c, s = transfer.ravel()
    mu = (a+s+1j*(c-b))/2
    nu = (a-s+1j*(c+b))/2
    ratio = -nu.conjugate()/mu.conjugate()
    # <0|M Q_j Q_i|0> / <0|M|0>, with j>=i.
    u = np.array([a + 1j*b for a, b in rows])
    later = np.tril((u.conj()[:, None]*u + ratio*u[:, None]*u)/2)
    covariance = later + np.tril(later, -1).T
    diff = np.diff(np.eye(n+1), axis=0)
    hessian = (diff.T@diff/d + np.diag(mass*np.array(curvature))).astype(complex)
    hessian[0,0] += 1j
    hessian[-1,-1] += 1j
    error = np.max(np.abs(covariance - np.linalg.inv(-1j*hessian)))
    print(n, error)
    assert error < 1e-10

if __name__ == "__main__":
    for n in (1, 2, 8, 32, 64):
        check(n)
