# Core computational optimization plan

Status: W1–W7 landed (2026-07-27). W7c landed safe GLRM gradient/proximal
refactor + `PalmSweepEnd` sweep context; full multi-block proximal reuse remains
blocked on per-oracle proof obligations. Lodged in mote as epic **Core spectral
and allocation optimization** (`bd-01KYK0E5JP9GADN2GS6J7BZ7PG`, closed). W7 parent:
`bd-01KYK0EMPR7H9A79Y2R0R20CT0` (closed).

## Goals

1. Reduce asymptotic spectral work on ordinary dense APIs (PCA, paired, SIMPLS,
   CPCA) without changing statistical definitions.
2. Eliminate the `copyData` → mutate → `matrixFromRowMajor` allocation tax that
   compounds across helpers.
3. Prefer Gale matvec/gemm/eigen selection over hand-rolled indexed loops.
4. Keep semantics and R/law tests green; floating-point differences must stay
   within existing tolerances.

## Non-goals

- Redesigning the typed semantic layer or package hierarchy.
- Premature micro-tuning of cold paths or outer `Either` validation.
- Shipping a new public estimator surface.

## Workstreams

### W1 — Paired identity / shared \(B^{-1/2}\) (P0)

**Problem.** PLSC/CCA/RRR pay full eigen (and sometimes SPD cert of `I`) to form
\(B^{-1/2}\) per side. Identity cometrics are still certified and inverted.

**Change.**
- Short-circuit identity normalization to `I`.
- Reuse SPD-cert eigen (or Cholesky) for ridge Grams instead of a second
  `inverseSquareRoot` eigen.
- Keep SPD rejection semantics.

**Verify.** Existing paired R fixtures and `PairedOperatorProblemSuite` /
`PairedLatentSuite`.

### W2 — Counted symmetric eigen + thin-Gram SVD (P0)

**Problem.** `DenseSolvers.symmetricEigen` always uses `EigenSelection.All`.
`GramSvdSolver` always forms \(X^\top X\) and fully diagonalizes it.

**Change.**
- Extend `SymmetricEigenSolver` with counted/leading selection.
- `GramSvdSolver`: choose thin Gram by `min(n,p)`; request only needed roots
  (with rank-cutoff padding as required).
- Keep full spectrum for SPD certs and callers that need \(\lambda_{\min}\).

**Verify.** PCA/SVD/CPCA suites, Gram SVD reconstruction laws, paired whitened
SVD.

### W3 — Allocation plumbing: `copyData` and builder sinks (P1)

**Problem.** `DMat.copyData` / `DVec.copyData` round-trip through Seq/Vector;
helpers rebuild via `matrixFromRowMajor` after Array mutation.

**Change.**
- Single contiguous primitive copy where Gale allows.
- Prefer writing into Gale builders; replace reinventions of diagonal-add and
  symmetrize with Gale primitives when available.
- Touch `GaleNumerics`, `MatrixOps`, then high-traffic call sites.

**Verify.** Core matrix/preprocess suites; no public API change.

### W4 — Dense matvec via Gale + SIMPLS buffer reuse (P1)

**Problem.** `DenseMatrixView` and SIMPLS hand-roll matvecs; SIMPLS allocates
per component and copies `S` on deflation.

**Change.**
- Route dense matvec through Gale `mulInto` / equivalent.
- Reuse SIMPLS work buffers; in-place scale/center; rank-1 update of `S`
  without full `copyData`.
- Prefer top-1 eigen for `computeQa` once W2 lands.

**Verify.** `PlsRegressionSuite` R `plsr` fixtures and score orthonormality.

### W5 — CPCA factored projectors + metric-root cache (P1)

**Problem.** Dense \(QQ^\top\) projectors and re-factoring of metrics; unused
blocks still materialize.

**Change.**
- Apply \(Q(Q^\top Z)\) without forming \(QQ^\top\).
- Cache metric roots from construction.
- Compute only requested blocks when diagnostics allow.

**Verify.** `CpcaOperatorProblemSuite` inertia partition and reconstruction.

### W6 — Column-affine closed-form inverse (P2) ✓

**Problem.** `predict`/`reconstruct`/`inverseContribution` ping-pong through
`inverseTransform` → materialize → densify for pure column affines; some paths
did it twice.

**Change (complete).**
- `FittedInvertiblePreprocessor.inverseTransformDense` / `inverseContributionDense`
  with closed-form overrides on `InvertibleColumnAffine` (full columns and column
  subsets).
- Wired through `FittedTransform.predict`, `FittedSynthesis` Original reconstruct,
  `Cpca.reconstruct`, and `EffectOperator.inverseContribution`.
- View + densify fallback retained for non-affine invertible prep.

**Verify.** Optimization-layer suites: `NumericsOptimizationSuite` (identity
short-circuit, counted eigen, Gram branches); `DecompositionSuite` row-Gram
(`n ≤ p`); `RowGeometrySuite` shared whitening + `selectRows`/`writeRows`;
`CanonicalEffectSuite` / `CanonicalSpectrumSuite` rank-limited leading roots;
`GpcaProblemSuite` cometric λ_min; `SparsePreprocessingSuite` dense-inverse;
`PalmConvergenceSuite` `PalmSweepEnd`; `CpcaOperatorProblemSuite` dense
reconstruct; plus existing law/R-parity suites.

### W7 — Follow-ons (P2) ✓

**W7a — Canonical/Fisher counted generalized eigen.**
- `fitSpectrum(k)` requests `Count(k)`; `fit()` requests
  `Count(min(n, rankEstimate(effect)))` for repeated leading roots without a full
  `n×n` pencil when the effect rank is smaller.
- LDA identity cometric uses analytic SPD certificate (paired W1 pattern).

**Verify.** `CanonicalEffectSuite` rank-limited `fit()` with repeated leading root;
existing canonical spectrum/oracle suites.

**W7b — GPCA Cholesky-backed cometric cert.**
- After Cholesky → inverse, λ_min certificate via `Count(1, SmallestAlgebraic)`
  instead of full `FormCertificates.spd` spectrum.

**Verify.** `GpcaProblemSuite` diagonal cometric λ_min oracle; R `genpca` parity.

**W7c — GLRM Θ cache / stationarity (safe subset).**
- `PalmSweepEnd` passes post-sweep update context into block stationarity probes.
- GLRM PALM: shared `proximalResidual` for update and stationarity; flat-array
  cache of row codes and decoder (Θ) during `entryGradient`; reused natural
  buffer across entries.
- Full reuse of update proximal maps at sweep end is **not** enabled: multi-block
  Gauss–Seidel stationarity must be measured at the final state.

**Verify.** `LowRankModelsFittedDifferentialSuite`, `PalmConvergenceSuite`,
`GeneralizedLowRankLifecycleSuite`.

**W7d — Block-Cholesky gather.**
- `RowGeometryOps.selectRows` / `writeRows` index elements directly; no full
  matrix `copyData`.

**Verify.** `RowGeometrySuite` block-Cholesky whitening and effect-operator suites.

**W7e — EffectModelFit shared whitening.**
- `EffectModelFit.fromTerms` whitens the design once and passes `Some(designW)`
  into per-term `EffectTermFit.fromDesign`.

**Verify.** `RowGeometrySuite` multi-term block-Cholesky projector equivalence test.

## Sequencing

```
W1 ──┐
     ├──► W4 (SIMPLS top-1 after counted eigen)
W2 ──┤
     └──► W5
W3 ──────► W4, W6 (allocation wins compound)
W7 after W1–W3 land and measure
```

## Acceptance

- `sbt compileAll testAll` green.
- Public surface unchanged unless explicitly decided (`tools/public-surface.sh`).
- Paired/PCA/PLS/CPCA law and R-parity suites remain within existing tolerances.
- Each workstream notes before/after cost narrative (what spectral or allocation
  work was removed), even without a formal benchmark harness.

## Risks

- Truncated eigen/SVD can change multiplicity handling; pad and document.
- Identity short-circuit must not skip intentional SPD evidence where the typed
  layer requires a certificate object — supply a trivial identity certificate.
- Allocation refactors must preserve MatrixView storage-policy semantics.
