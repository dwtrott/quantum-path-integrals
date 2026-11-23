# Working rules for AI collaborators

These rules apply to Nova and Sparks. David coordinates task assignment and acceptance.

## Claims

1. Distinguish three kinds of support in every report: numerical evidence, conventional proof, and Lean proof that has actually compiled.
2. Never report a compiler run that did not occur. A Lean result counts as verified only when CI passes on the pull request that contains it.
3. No accepted theorem depends on `sorry`, `admit`, or an unapproved axiom. Run `#print axioms` on each accepted theorem and check the statement against the mathematical claim.
4. Floating point agreement is evidence. Numerical certification requires interval arithmetic or another validated method.
5. Do not claim novelty without a literature check, and not because a formalization is absent.
6. State for every result which object it concerns: the regulated finite-dimensional integral, the finite-dimensional integral with the regulator removed, or the continuum limit.

## Workflow

1. Work on a branch named `nova/...` or `sparks/...` and merge through a pull request reviewed by David.
2. Keep the toolchain and dependency pins committed: `lean-toolchain`, `lakefile.toml`, `lake-manifest.json`.
3. Record parameters, commands, library versions and the commit hash with each numerical result.
4. Keep failed attempts, labeled as such, under `research/counterexamples/` or in the relevant report.
5. Record each proposed insight, proof or counterexample in `research/ledger.md` with date, contributor, reviewer and verification status.
6. Do not impersonate the other collaborator. Do not commit credentials or unrelated private material.

## Layout

- `QuantumPathIntegrals/`: Lean source
- `proofs/`: conventional, human-readable arguments
- `experiments/`: scripts and reproducible outputs
- `research/`: correspondence, literature notes, conjectures, counterexamples, ledger
- `reports/`: benchmark and proof-review reports
