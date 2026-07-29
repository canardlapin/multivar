# TaskComponents methodology plan

Status: merge gate accepted (2026-07-27). Phase 0–3 landed.
Ownership: `family.paired` owns the estimand; `family.spectral` owns
metric-coordinate execution; `family.cpca` keeps crossed G×H blocks.
No `family.task` for v1.

Related: [`multivar-package-hierarchy.md`](multivar-package-hierarchy.md),
[`paired-latent-backbone.md`](paired-latent-backbone.md),
[`multivar-public-api-surface.md`](multivar-public-api-surface.md).

## One-liner

```text
TaskComponents owns a new estimand, not a new SVD stack.
```

The estimand is the low-rank structure and quadratic gain of a least-squares
task after resolving it into a scientifically specified focus and its
orthogonal remainder.

Intellectual sources (documentation only; not public API switches):

- De la Torre — least-squares reduced-rank task.
- Allen–Grosenick–Taylor — generalized metric decomposition of the canonical
  operator.
- Takane–Hunter — weight / structure / covariance / regression resolution
  semantics.

## Settled ownership

```text
family.spectral
  MetricRoots (existing)
  MetricCoordinates (extract)
  OrthonormalSubspace (extract)
  LocalSvd / exact GMD execution
          │
          ▼
family.paired
  CanonicalTask
  TaskResolution
  DualSpectrum
  RankBudget
  TaskComponents
  TaskValue / TaskComponentBlock
  TaskComponentsFit
  TaskDecompositionFit (Phase 3)

family.cpca          (unchanged)
  Cpca, GxH / G0xH / GxH0 / G0xH0
```

| Package | Owns | Does not own |
| --- | --- | --- |
| `family.spectral` | metric roots, whitened coordinates, range extraction, local SVD/GMD | task resolution vocabulary, value accounting |
| `family.paired` | resolved LS-task estimand, resolution grammar, predictive fit, gain accounting | crossed table blocks; a third metric-SVD stack |
| `family.cpca` | two-sided row×feature block CPCA | Takane weight/structure on K |
| `analysis` | curated public names | implementation types |

Reverse imports into `family.task` are moot: there is no `family.task`.
`family.spectral` and `family.cpca` must not import TaskComponents types.
`family.paired` may depend on spectral substrate only.

## Estimand axes

```text
Task × Geometry × Resolution
```

- **Task** — what is approximated or predicted (De la Torre).
- **Geometry** — how error and component magnitude are measured (GMD).
- **Resolution** — which scientific part of the task is separated from its
  orthogonal remainder (Takane–Hunter).

## Canonical operator

For the weighted prediction task \(Y \approx XB\) with observation metric \(W\),
target metric \(R\), and (when allowed) source ridge \(\rho\):

\[
G_\rho = X^\top W X + \rho I,\qquad
C = X^\top W Y,\qquad
K = G_\rho^{-1/2} C R^{1/2}.
\]

Phase 1 uses certified SPD inverse square roots (current paired policy), not
Moore–Penrose \(G^+\). DualSpectrum policies apply only to the singular
spectrum of \(K\), never as a silent substitute for Gram normalization.

For a decomposition self-task (Phase 3):

\[
K = Q^{1/2} X R^{1/2}.
\]

This \(K\) is already latent in `PairedOperatorProblem.solve` as

```text
whitened = sourceInverseHalf * scaledCross * targetInverseHalf
```

### Metric vs normalization (must not conflate)

Public vocabulary: \(R\) is the **target metric** in the loss
\(\operatorname{tr}[W(Y-XB)R(Y-XB)^\top]\).

The paired substrate may continue to apply a **target normalization** \(T\) via
\(T^{-1/2}\). Representing loss metric \(R\) requires \(T = R^{-1}\). For CCA,
\(T = \Sigma_{YY}\) and \(R = \Sigma_{YY}^{-1}\), so
\(T^{-1/2} = R^{1/2}\).

### Mathematical gate: ridge × nonidentity \(R\)

Existing Frobenius ridge \(\rho\|B\|_F^2\) with \(R \neq I\) yields the
Sylvester equation

\[
X^\top W X\, B R + \rho B = X^\top W Y\, R,
\]

which is **not** the one-sided whitening problem that produces \(K\) above.
Target-compatible ridge \(\rho\,\operatorname{tr}(B R B^\top)\) is a later
named contract.

**Phase 1 Options must allow only:**

1. identity target metric + existing `RegressionRegularization` ridge; or
2. arbitrary target metric with \(\rho = 0\) (OLS).

Reject other combinations with a typed error that names the Sylvester gap.

## Resolution grammar

```scala
TaskResolution.Whole
TaskResolution.Target.weights(design)
TaskResolution.Target.structure(design, dual = DualSpectrum.moorePenrose)
TaskResolution.Source.covariance(design)   // Phase 2
TaskResolution.Source.regression(design)   // Phase 2
```

Cases are private; constructors are grammatical. No `useDual` booleans and no
public A–D letter codes. Letters belong only in mathematical documentation.

`DualSpectrum`:

- `MoorePenrose(rankTolerance)`
- `Truncated(components, rankTolerance)`
- `Ridge(lambda)`

These are distinct estimands for structure resolution, not interchangeable
numerical backends.

### Axis meaning (document from Phase 1)

| Task | Source axis of \(K\) | Target axis of \(K\) |
| --- | --- | --- |
| Decompose \(X\) | observations (rows) | features |
| Predict \(X\to Y\) | predictor features | response features |

Phase 2 `Source.*` designs therefore have one row per predictor feature in
prediction mode. Shape errors must name the axis, not only “matrix shape
mismatch.”

## RankBudget

```scala
sealed trait RankBudget
object RankBudget:
  final case class Focus(components: Int) extends RankBudget
  final case class Total(components: Int) extends RankBudget
  final case class Split(focus: Int, remainder: Int) extends RankBudget
```

| Plan | Semantics |
| --- | --- |
| `Focus(k)` | \(k_f = k\), \(k_r = 0\) |
| `Split(k_f, k_r)` | explicit branch ranks |
| `Total(k)` | \(k_f + k_r \le k\), largest local \(d_{bj}^2\) across branches |

Ordinary shorthand `components = k` means `RankBudget.Focus(k)`.
`Total` must report whether the selection boundary cuts a tied spectral
cluster.

## Public Phase 1 surface

```scala
import multivar.analysis.*

val fit =
  TaskComponents.fit(
    predictors = x,
    responses = y,
    resolution =
      TaskResolution.Target.structure(
        design = h,
        dual = DualSpectrum.Ridge(lambda = 1e-4)
      ),
    rank = RankBudget.Focus(components = 4)
  )
```

Curated `analysis` exports:

```scala
TaskComponents
TaskComponentsFit
TaskResolution
DualSpectrum
RankBudget
```

`TaskValue` / `TaskComponentBlock` remain public return types under
`family.paired` but need not be named façade exports. No public
`CanonicalTask`, projector algebra, or `TaskAlgebra`.

Result vocabulary (loss language, not “scientific worth”):

```text
fit.predict / coefficients / intercept
fit.focusContribution / remainderContribution
fit.focus / remainder   (scores, patterns, singular values, per-part gains)
fit.value.baselineLoss / fullLoss / fittedLoss
fit.value.availableGain / retainedGain / unretainedGain
```

Exact identities:

\[
\text{availableGain}
=
\text{focus.availableGain}
+
\text{remainder.availableGain},
\]

\[
\text{fittedLoss}
=
\text{baselineLoss}
-
\text{retainedGain},
\]

\[
\widehat Y
=
\text{baseline}
+
\widehat Y_f
+
\widehat Y_r.
\]

`Whole` yields a lawful empty remainder (`rank == 0`, `availableGain == 0`).
Ordinary documentation must not lead with Whole as an alternative spelling of
RRR.

## Chooser contrast (site-docs)

| Object | External structure | Method |
| --- | --- | --- |
| Data table | none | Pca / Gpca |
| Data table | crossed row × feature blocks | Cpca |
| Directed prediction task | none | ReducedRankRegression |
| Canonical LS task operator \(K\) | focus ⊕ orthogonal remainder | TaskComponents |

Documentation sentence (required):

> Cpca partitions table inertia into crossed G×H blocks. TaskComponents
> resolves a least-squares task operator according to weight, structure,
> covariance, or regression semantics. Neither is a replacement for the other.

## Non-goals (v1)

- Arbitrary Hadamard entry weights, missing-data loss, sparse/nonnegative
  factors, resolution chains, kernel builders.
- Automatic selection between weight and structure semantics.
- Replacing or wrapping public `Gpca` / `Cpca` / `ReducedRankRegression`.
- Silent substitution of DualSpectrum ridge for ill-conditioned Moore–Penrose.
- Pseudoinverse Gram normalization (`NormalizationSupport.EffectiveRange`)
  until quotient-space behavior is specified for paired, GPCA, decoding, and
  gain accounting.

## Phases

### Phase 0 — substrate extraction (no public API)

Extract or consolidate under `family.spectral` / `capability`:

- `MetricCoordinates` with `table` and `cross` constructors.
- `OrthonormalSubspace` (`Z ↦ E(EᵀZ)`, no materialised \(EE^\top\)).
- `LocalSvd` (or thin wrapper over existing `SvdSolver`).
- `capability.PairedCoordinateMap` — shared raw coefficient/intercept decode
  used by RRR and TaskComponents.

Refactor `PairedOperatorProblem.solve` onto `MetricCoordinates.cross` without
changing PLSC / CCA / RRR outputs. Prefer opportunistic reuse for GPCA/CPCA
whitening paths; parity suites must stay green.

**Exit:** `sbt compileAll testAll`; public surface unchanged.

### Phase 1 — predictive target resolution

Add under `family.paired`:

- `TaskResolution` (Whole, Target.weights, Target.structure)
- `DualSpectrum`
- `RankBudget`
- `CanonicalTask`, `TaskResolver`, `SpectralTransport` (package-private)
- `TaskComponents`, `TaskComponentsFit`, `TaskValue`, `TaskComponentBlock`

Public support: Focus / Total / Split; R×ridge Options gate; axis metadata on
fits.

**Laws:**

- RRR equivalence under Whole (identical prep, ridge, SPD, solver).
- Focus/remainder orthogonality of \(K\).
- Available-gain additivity; fitted-loss identity.
- Contribution reconstruction identity.
- Design equivariance \(H \mapsto HA\) for nonsingular \(A\).
- DualSpectrum provenance: Moore–Penrose / Truncated / Ridge labelled distinctly.

**Docs:** `site-docs/analyses/task-components.md`; chooser row; hierarchy note.

**Exit:** façade export + `tools/public-surface.sh --update` reviewed;
`docsCheck` green.

### Phase 2 — source resolution

Add `Source.covariance` and `Source.regression` with axis-aware shape errors
and comparison tests. Still predictive-only.

### Phase 3 — resolved decomposition with inductive feature transfer

`TaskComponents.decompose` with nontrivial resolution required on the ordinary
path. Unresolved metric decomposition remains `Gpca`. Lowers through shared
spectral coordinates and transports left resolutions to feature-side maps
\(F = Z^{+}EZ\) so new rows support:

- `transform` / `reconstruct`
- `focusContribution` / `remainderContribution` / `residual`
- `TaskResolutionTransfer` (`FeatureSpecified` vs `TrainingSourceInduced`)

**Laws:** training-score reproduction; training reconstruct identity; new-data
contribution assembly; \(ZF_b \approx Z_b\); approximate transfer idempotence;
target design equivariance on new-data scores; GMD/PCA Whole only as an
internal law.

**Exit:** `TaskDecompositionFit` inductive API exported; public Whole rejected;
docs/chooser updated; `compileAll testAll docsCheck` green.

## Merge checklist

- [x] Hierarchy plan amended (this document + package hierarchy bullet).
- [x] Phase 0 substrate extracted; PLSC/CCA/RRR/(GPCA/CPCA) parity green.
- [x] Phase 1 Options reject nonidentity \(R\) with Frobenius ridge.
- [x] RRR Whole law + gain/contribution identities green.
- [x] Chooser documents Cpca ≠ TaskComponents.
- [x] Public surface review for new analysis exports.
- [x] No `family.task` package created.
- [x] Phase 3 `decompose` rejects Whole; GMD Whole law is internal-only.
- [x] Phase 3 inductive `transform` / `reconstruct` via \(F = Z^{+}EZ\).
- [x] Follow-on: biorthogonal `FittedBidirectionalTransform`, branch transform, `validationLoss`.
- [x] Follow-on: `resolve` / `unsupportedFraction` (\(F_0 = I - P_{Z^\top}\)).

## File layout (target)

```text
modules/core/shared/src/main/scala/multivar/
  family/spectral/
    MetricCoordinates.scala          # Phase 0
    OrthonormalSubspace.scala        # Phase 0
  capability/
    PairedCoordinateMap.scala        # Phase 0
  family/paired/
    TaskResolution.scala             # Phase 1
    DualSpectrum.scala               # Phase 1 (or co-located)
    RankBudget.scala                 # Phase 1
    TaskComponents.scala             # Phase 1
    TaskComponentsFit.scala          # Phase 1
    TaskValue.scala                  # Phase 1
    internal/                        # optional package-private grouping
      CanonicalTask.scala
      TaskResolver.scala
      SpectralTransport.scala

modules/core/shared/src/test/scala/multivar/family/paired/
  TaskComponentsSuite.scala
  TaskValueLawsSuite.scala
  TaskReductionLawsSuite.scala

site-docs/analyses/task-components.md
```

## Relationship to existing methods

Familiar direct methods remain. TaskComponents is compositional for questions
such as:

- predictive components of \(X\to Y\) resolved by response-domain groups;
- weight-constrained vs structure-constrained interpretations of the same \(H\);
- partition of task gain into a prespecified feature system and its complement.

Some existing methods may later reuse CanonicalTask / MetricCoordinates
internally. That is an implementation benefit, not a reason to abstract their
public APIs.
