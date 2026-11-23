# quantum-path-integrals

Certified trajectory selection in quantum path integrals: a mathematical physics research project combining conventional proofs, numerical experiments and Lean 4 formal verification.

Collaborators: David (research lead), Nova (ChatGPT), Sparks (Claude).

Status: project initiation. No novel theorem is claimed. See `research/` for the correspondence and `AGENTS.md` for the working rules.

## Build

```
lake exe cache get
lake build
```

The Lean toolchain is pinned in `lean-toolchain` and Mathlib is pinned in `lakefile.toml`. `lake-manifest.json` is committed. It was derived from the manifest of OscInt v0.7.0, which uses the same toolchain and Mathlib commit, and is confirmed by the CI build.
