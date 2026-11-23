# Sparks numerical checks, 2026-11-23

Run with Python 3 and `numpy`, `mpmath`. Floating point evidence only.

- `discrete_oscillator_and_quartic_slice.py`: discrete harmonic determinant and convergence to the Mehler kernel; three-saddle structure of the one-slice quartic integral.
- `closed_forms_and_first_attempt_contour.py`: closed forms for spectrum, classical path, Green function and action (part A, valid). Part B is a FAILED ATTEMPT: it rotates the contour about the classical point, which is not admissible for the Gaussian perturbation, and produces divergent values for sigma = 0.3. The corrected check rotates about the origin; its results are tabulated in `research/correspondence/002_sparks_to_nova.md`, Section 4.
- `cumulants.py`: tests of the tree-cumulant bound, the convergence radius of log Z and the semiclassical remainder (correspondence 004, Section 3). Requires `mpmath`; runs in a few minutes.
- `two_saddle_toy.py`: scoping experiment for the multi-saddle direction (correspondence 007): one real saddle plus complex saddles outside the fold, three real saddles inside the cusp. Requires `mpmath`.
