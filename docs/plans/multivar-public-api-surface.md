# Multivar public API surface plan

Status: Phases 0–5 landed for local/CI. Phase 4 topology uses a Gale Maven
coordinate + `publishLocal` pin, consumer smoke, and MiMa scaffolding; a public
Maven Central release still waits on Gale's own Central publication. Phase 5
added SIMPLS PLS regression behind the ordinary façade as the first post-0.1
architectural vertical. Phase 6 onward remains deferred feature work. This
document responds to an
external API-surface review of the 0.1 line. It records which review claims
were verified against the code, where this plan deliberately departs from the
review, and the order in which the work should be done.

Scope note. This is a *surface* plan. It does not revise the mathematical
contract in [`multivar-mathematical-contract.md`](multivar-mathematical-contract.md)
or the review obligations in [`multivar-external-review.md`](multivar-external-review.md),
both of which concern what multivar may claim. This document concerns what an
ordinary caller is allowed to see and forced to preserve. It amends the import
policy in [`multivar-package-hierarchy.md`](multivar-package-hierarchy.md).

## Diagnosis

The reviewer's verdict is correct and is accepted: the core architecture is
sound and should not be simplified; the public boundary around it is not
selective enough to publish. Two consequences drive everything below.

1. Every implementation type reachable from a fitted result becomes a
   compatibility obligation the moment 0.1 is tagged. Today `SvdResult`,
   `PairedOperatorFit`, `GpcaOperatorFit`, `LdaOperatorFit`, `NystromOperatorFit`,
   `FittedFrameTransform`, `KernelEigenArtifact`, and `NystromState` are all
   reachable from ordinary fits.
2. Every ordinary fit is a `case class`, so callers get a public `apply`,
   `copy`, and `unapply`. The relationships established during fitting are
   therefore not enforced by the type: a caller can pair one `SvdResult` with
   an unrelated `FittedFrameTransform`, or replace one side of a paired fit.
   The deep core is typed; the artifact handed to users is forgeable.

The single highest-value change is therefore to make fitted results opaque
before release, because that is what determines how much of the current
implementation must be preserved indefinitely.

## Verification of review claims

Each claim was checked against the code at `81494d8`.

| Claim | Verdict | Evidence |
| --- | --- | --- |
| Ordinary fits are public case classes leaking internals | Confirmed | `PcaFit(result: SvdResult, transform: FittedFrameTransform)` at `family/spectral/Decompositions.scala:47`; `NystromFit` has 13 fields at `family/kernel/Kernel.scala:292` |
| Existential `SemanticSpace` parameters appear in public signatures | Confirmed | `GpcaFit.result`, `PlscFit.operator`, `CcaFit.operator`, `ReducedRankRegressionFit.operator`, `LdaFit.result` all take `[? <: SemanticSpace, ...]` |
| RRR `coefficients` are not in original coordinates | Confirmed | coefficient is built from working-space factors at `family/paired/PairedDecompositions.scala:346-351` and stored unchanged by `FittedCoefficientTransform.from`; no `B_raw = D_x B_w D_y^{-1}` path exists |
| Documentation promises original coordinates | Confirmed | `site-docs/analyses/paired.md:78` |
| No intercept is exposed for RRR | Confirmed | no `intercept` member anywhere; the shift lives implicitly inside the fitted response preprocessor |
| GPCA applies ordinary centering before non-ordinary geometry | Confirmed | 4-arg `Gpca.fit` defaults to `PreprocessSpec.Center` at `family/spectral/Gpca.scala:42`, applied before the row metric reaches `DynamicGpcaProblem` |
| The semantic layer already has a richer centering concept | Confirmed | `CenteringPlan` with `NoCentering`, `CenterByMeasure`, `CenterOrthogonally`, `AlreadyCenteredBy` at `core/SemanticDiagram.scala:381` |
| CPCA does no preprocessing | Confirmed | `Cpca.fit` passes the raw input to `CpcaOperatorProblem.fromMatrices`; `CpcaEstimatorSpec` has no preprocessing field |
| `Lda` is a projection, not a classifier | Confirmed | `LdaFit` exposes only `scores`, `weights`, `criterionValues`, `project`; no classes, priors, centroids, decision function, or `predict` |
| `ridge` is a trace-scaled fraction | Confirmed | `TraceRidgeFraction` wrapped at `family/canonical/Lda.scala:596`, scaled by `trace(W)/dim` at line 493 |
| Nyström silently alters the landmark request | Confirmed | `indices.toVector.distinct.sorted` runs before validation at `family/kernel/Kernel.scala:159`, while `MultivarError.DuplicateIndex` already exists and is used by `IndexSet.from` |
| The landmark kernel is symmetrized unconditionally | Confirmed | `MatrixOps.symmetrize(kMmRaw)` at `family/kernel/Kernel.scala:429`; no asymmetry is measured; `Kernel` is an open trait |
| `PreprocessSpec` naming and invertibility are imprecise | Confirmed | one `FittedPreprocessor` trait carries both directions; a zero scale fails only at `inverseTransform` (`MultivarError.NonInvertibleValue`), so RRR can fit and then fail to predict |
| No curated façade exists | Confirmed | no `multivar.analysis`, `multivar.advanced`, or `multivar.syntax` package; no `orThrow` |
| Adding an estimator is non-local | Confirmed, narrower than stated | roughly 8–10 exhaustive `match` sites in `workflow/Plans.scala`, plus `PairedProgramKind` sites; the IR is *not* implicated, because operator plans deliberately do not serialize the estimator enum (`ir/.../OperatorPlanIrSuite.scala:24`) |
| Gale is a pinned Git source dependency | Confirmed | `ProjectRef` on `github.com/canardlapin/gale.git#d55fe2f...` at `build.sbt:37-43` |
| No MiMa, no consumer smoke project | Confirmed | `project/plugins.sbt` carries only scalajs, crossproject, and `sbt-typelevel-site` |

Two review statements need correction.

**The data for `explainedVarianceRatio` are not already available.** The review
asserts the omitted PCA quantities can be computed from what is retained. That
holds for `explainedVariance`, which is `s_i^2 / (n - 1)`. It does not hold for
the *ratio*, which needs the total variance of the preprocessed training data.
`TransformDiagnostics.spectrum` stores only the retained singular values and is
required to have length `effectiveComponents` (`capability/FittedTransform.scala:9-17`),
so for a truncated fit the denominator is unrecoverable. Total sum of squares
must be captured during fitting. Likewise `center` and `scale` are not on the
`FittedPreprocessor` trait, only on the concrete `FittedColumnAffine`, so
exposing `pca.center` needs a new narrow accessor rather than a delegation.

**A wildcard export from a package is illegal, but a curated named export is
not.** The hierarchy doc rejects root aliases partly on the grounds that a
package cannot be wildcard-exported. That restriction is real but narrow: Scala
3 forbids a *wildcard or given* selector whose qualifier is a package, and
permits top-level export clauses with named selectors. The proposed façade uses
named selectors exclusively, so it compiles as written and needs no hand-
maintained alias per member. The policy amendment below is therefore a genuine
editorial decision, not a workaround.

## Departures from the review's plan

**Fix semantics and redesign results per family, not as two global passes.**
The review sequences "redesign every fitted result" ahead of "fix the semantic
contract issues". Taken literally that touches each result type twice, because
the RRR fix changes which members exist at all: `coefficients` changes meaning,
`intercept` appears, and `fullCoefficient` is renamed. The plan below does one
family at a time, landing its opaque result and its semantic correction in the
same change.

**Preprocessing precision comes first, not sixth.** Deriving raw-coordinate RRR
coefficients requires `D_y^{-1}` at fit time, which requires the invertible
preprocessor distinction. It is a prerequisite for the paired family, not a
later tidying step.

**`GpcaCentering.Auto` refuses to guess rather than dispatching silently.**
The review proposes an `Auto` case that selects ordinary centering, a weighted
mean, or metric-orthogonal removal according to the row metric. That
reintroduces exactly the hidden estimand-changing policy choice the same
section objects to, just with more branches. Instead, `Auto` means ordinary
centering when the row metric is the identity and a typed error otherwise. The
convenient path stays convenient for the common case; a nonuniform row metric
forces the caller to name the centering they want.

**The name `Lda` is not reused.** The review suggests renaming the projection
and later adding `LdaClassifier`. If `Lda` were vacated and later reintroduced
with classifier semantics, old code and old documentation would silently change
meaning. `Lda` is retired outright: `FisherDiscriminant` for the projection,
`LdaClassifier` for the classifier.

**CPCA keeps one name.** Rather than splitting `Cpca.fit` and
`ConstrainedSvd.fit`, `Cpca.fit` gains a `preprocessing` parameter defaulting
to `PreprocessSpec.Center`. The as-supplied behaviour remains reachable as
`preprocessing = PreprocessSpec.Pass`, which is what a second name would have
meant anyway, without adding a second estimand to document.

**No deprecation shims.** The README already states that APIs and result types
may change without a migration period, and nothing is published. Renames
should be clean. This is the last cheap opportunity to take them.

## Target ordinary surface

Three rings, enforced by the package surface rather than only described in
prose.

| Ring | Package | Audience | Contents |
| --- | --- | --- | --- |
| 1 | `multivar.analysis`, `multivar.syntax.unsafe` | most analysts | dense estimators, simple configuration, opaque fitted results |
| 2 | `multivar.core`, `multivar.family.*`, `multivar.workflow` | advanced Scala users | `MatrixView`, schemas, metrics, plans, storage policies, checked request types |
| 3 | `multivar.advanced` and the layers it opens | library authors, mathematical extensions | operators, typed frames, certificates, provenance, solver receipts, lowering |

Ring 3 is reached only through an explicit import. Extension methods in
`multivar.advanced` mean `pca.typedFrame` does not appear in autocomplete until
the caller asks for it.

The intended first page:

```scala
import multivar.analysis.*

val pca = Pca.fit(x, components = 5, preprocessing = Preprocessing.Center)
val fit = pca.getOrElse(sys.error("fit failed"))
println(fit.explainedVarianceRatio)
```

`fit` continues to return `Either[MultivarError, PcaFit]`. `multivar.syntax.unsafe`
adds `orThrow` for worksheets and demonstrations; exceptions never become the
primary interface.

### Verb convention

| Verb | Meaning |
| --- | --- |
| `transform` | map observations to the fitted latent representation |
| `inverseTransform` | map latent coordinates back to feature space where a fitted reconstruction contract exists |
| `predict` | estimate directed responses |
| `reconstruct` | reconstruct supplied or training data |
| `encode` | infer latent variables where encoding is itself an optimization problem |
| `project` | retained on the typed/advanced layer as a mathematical alias |

Ordinary results therefore use `transform`, `transformX`/`transformY`,
`predict`, `reconstruct`, `encode`. `project` leaves the ordinary surface and
stays on `FittedFrameTransform` and the typed frames, where it names an
operator action rather than a user task.

## Work plan

### Phase 0 — Instrumentation (do first)

Nothing else is measurable until the surface is visible. The review places
public-surface tests sixth; the measuring instrument belongs first, because it
is what tells us whether the later phases actually worked.

- **0.1 Signature snapshot. Landed.** `tools/public-surface.sh` dumps the public
  members of the ordinary surface with `javap -public`, normalized and sorted,
  into `tools/public-surface/ordinary-api.txt`, and fails when the committed
  file is stale. The surface it covers is listed explicitly in
  `tools/public-surface/surface-classes.txt`, so adding a type to the ordinary
  API is a deliberate edit rather than a side effect. CI runs `--check` after
  compilation. This is the pre-release substitute for MiMa, and it makes the
  leakage described above visible as a diff. Because it records the surface as
  it is rather than asserting a target, it is the one instrument that can be
  built before any of the redesign; the first snapshot is a baseline of the
  problem, not of the goal.
- **0.2 Documentation correctness fix. Landed.**
  `site-docs/analyses/paired.md` promised that RRR `coefficients` are in
  original predictor and response coordinates. That is false under any scaling
  preprocessing, and the default is `Center`. The text now states the working
  coordinates, the centering-only special case in which the slopes coincide,
  the implicit intercept, and the instruction to use `predict` for predictions
  in original units.

The negative compile-time surface guard — asserting that `SemanticSpace`, `Op`,
`OperatorProgramFit`, `SvdResult`, `FittedFrameTransform`, certificates, and
path-dependent core types are unreachable from a file importing only
`multivar.analysis.*` — cannot be written before the façade and the opaque
results exist, since today every one of those assertions would fail. It is
therefore scheduled in Phase 3 alongside the façade rather than here. The
snapshot covers the intervening period: during Phase 2 the leaked types leave
the golden file one family at a time, which is a reviewable record of progress.

The baseline snapshot is worth reading once as a statement of the problem. It
shows `apply`, `copy`, and `unapply` on every fitted result; `PcaFit.result`
returning `SvdResult`; `? extends SemanticSpace` wildcards on the GPCA, paired,
and LDA fits; and thirteen constructor parameters on `NystromFit`, among them
`NystromState`, `KernelEigenArtifact`, and `NystromOperatorFit`.

### Phase 1 — Preprocessing foundation — **landed**

Prerequisite for the paired family.

- **1.1 Split invertibility into the type.** *Done.* `FittedPreprocessor` keeps
  `transform` and `restrict`; `FittedInvertiblePreprocessor` adds
  `inverseTransform`. Response-side preprocessors in regression signatures
  require the invertible type: `FittedCoefficientTransform` now stores and
  demands one, and the RRR fit path calls `fitInvertible` on the response spec,
  so a zero response scale is rejected at fit time instead of surfacing as a
  failed `predict`.
- **1.2 Rename and refine the spec.** *Done.* `Scale(weights)` is
  `MultiplyColumns(weights)`, and `PreprocessSpec.scale` is
  `PreprocessSpec.multiplyColumns`. `Standardize(convention)` carries a
  `VarianceConvention` defaulting to `Sample`, and `ColumnStats` gained
  `standardDeviations(convention)` so `Population` is a real alternative rather
  than a label. The scikit-learn-compatible treatment of degenerate columns and
  its test are unchanged.
- **1.3 Narrow affine summary.** *Done.* `ColumnAffineSummary(center, scale)`
  reports the fit in the familiar `(x - center) / scale` form, available from
  `InvertibleColumnAffine.summary`. It is not yet reachable from any fitted
  result; Phase 2 decides how each family exposes it.

Two departures from the plan as written, both toward a more static guarantee:

- `PreprocessSpec.fit` still returns `FittedPreprocessor`. Returning "the
  invertible type when the scale is invertible" can only vary the runtime class,
  not the static one, so callers would have to pattern-match or downcast to use
  it — and `family/kernel/Kernel.scala` already matches on the concrete
  `FittedColumnAffine` to derive its centering, which such a change would
  silently break. Invertibility is instead requested explicitly:
  `FittedPreprocessor.requireInvertible` returns the proof or a typed error, and
  `PreprocessSpec.fitInvertible` is the one-step form used by fit paths.
- The proof is a value, not a check. `InvertibleColumnAffine` holds the inverse
  affine computed at construction, so `inverseTransform` cannot fail on a zero
  scale — undoing `x * scale + shift` is applying `x * (1/scale) +
  (-shift/scale)`, which is another affine and needs no separate machinery.

Callers that legitimately invert a preprocessor they did not fit —
`FittedSynthesis` reconstructing in original coordinates, and `RowGeometry`
computing an effect contribution — now ask for the proof and propagate the typed
error. `BlockwisePreprocessor` implements `requireInvertible` by delegating to
its global affine.

Files changed: `core/Preprocessing.scala`, `core/MatrixView.scala`,
`core/RowGeometry.scala`, `capability/FittedTransform.scala`,
`capability/FittedSynthesis.scala`, `family/multiblock/Multiblock.scala`,
`family/paired/PairedDecompositions.scala`,
`validation/RecoveryValidation.scala`. Tests:
`core/SparsePreprocessingSuite.scala`, `family/multiblock/MultiblockSuite.scala`,
plus mechanical `Standardize()` updates in the workflow and validation suites.
`compileAll testAll` passes on JVM and Scala.js; `docsCheck` passes.

### Phase 2 — Opaque results and semantic fixes, family by family

Each item lands one family's private-constructor result, its statistical
accessors, its verb normalization, its advanced escape hatch, and its semantic
correction together.

- **2.1 Spectral: `Pca`, `Svd`.** *Landed.* `PcaFit` and `SvdFit` are
  `final class`es with `private[multivar]` constructors. Ordinary accessors are
  `scores`, `loadings`, `components`, `singularValues`, variance/inertia
  shares, `center`/`scale`, component counts, `transform`, `inverseTransform`,
  and `reconstruct`. Fitting captures the preprocessed total sum of squares so
  truncated ratios keep a stable denominator. PCA reports
  `explainedVariance` / `explainedVarianceRatio`; SVD reports `inertia` /
  `inertiaRatio`. Shared arithmetic lives in package-private `SpectralCore`
  that never appears in a public constructor signature (Scaladoc crashes on
  that pattern). `diagnostics` remains on the typed frame and is opened through
  `multivar.advanced.typedFrame` rather than on the ordinary fit.
- **2.2 Paired: `Plsc`, `Cca`, `ReducedRankRegression`.** *Landed.* Opaque
  results with `transformX`/`transformY`, `correlations`/`covariances`, and the
  paired weights and loadings. RRR exposes `workingCoefficients`,
  `coefficients` (raw predictor to raw response), and `intercept`.
  `unconstrainedWorkingCoefficients` lives on `multivar.advanced`.
- **2.3 GPCA centering.** *Landed.* `GpcaCentering` replaces a generic
  `PreprocessSpec` default. `Auto` resolves to ordinary centering under an
  identity row metric and returns a typed error otherwise. Opaque `GpcaFit`
  stores the training feature centre.
- **2.4 Canonical: retire `Lda`.** *Landed.* `FisherDiscriminant.fit` accepts
  generic ordered labels and takes an explicit `WithinScatterPolicy`
  (`FixedTraceScaledRidge` / `RequirePositiveDefinite`). Opaque
  `FisherDiscriminantFit` uses `transform`. Soft incidence and nuisance remain
  on the operator path; `LdaClassifier` is deferred.
- **2.5 CPCA preprocessing.** *Landed.* `Cpca.fit` and `CpcaEstimatorSpec`
  default to `PreprocessSpec.Center`. Reconstruction distinguishes
  `reconstructWorking` (metric / preprocessed space), `reconstruct` (original
  features after inverse preprocessing), and block
  `reconstructWhitened` / `reconstructMetric`.
- **2.6 Nyström request fidelity.** *Landed.* `LandmarkSet.from` rejects
  duplicates; `canonicalize` is explicit. `KernelSymmetryPolicy` defaults to
  measuring asymmetry and rejecting material violations. Opaque `NystromFit`
  exposes scores, eigenvalues, component counts, `transform`, and diagnostics.
- **2.7 GLRM.** *Reviewed.* Already private-constructor and lifecycle-wrapped;
  `CanEncode` / `FittedLatentEncoder.encode` remain the advanced encode path.
  No ordinary façade export.

Files: `family/spectral/Decompositions.scala`, `family/spectral/Gpca.scala`,
`family/paired/PairedDecompositions.scala`, `family/paired/PairedOperatorProblem.scala`,
`family/canonical/Lda.scala`, `family/cpca/Cpca.scala`, `family/kernel/Kernel.scala`,
`workflow/Plans.scala`. Tests: the corresponding suites plus
`family/spectral/DecompositionSuite.scala`, `family/paired/PairedLatentSuite.scala`,
`family/canonical/LdaSuite.scala`, `family/cpca/CpcaOperatorProblemSuite.scala`,
`family/kernel/KernelSuite.scala`, `workflow/PlanSuite.scala`.

Size: large. Phase 2 is landed.

### Phase 3 — Façade and documentation — **landed**

- **3.0 Surface guard suite.** *Done.* `AnalysisSurfaceSuite` uses
  `typeCheckErrors` to assert that `import multivar.analysis.*` resolves the
  dense estimators and preprocessing vocabulary, does not bring
  `SemanticSpace`, `Op`, `OperatorProgramFit`, `SvdResult`,
  `FittedFrameTransform`, or `VariationalFrameCertificate` into scope, and does
  not resurrect `multivar.Pca`.
- **3.1 `multivar.analysis`.** *Done.* Named-selector exports of the dense
  estimators and fits, `PreprocessSpec`, `VarianceConvention`, and
  `MultivarError`. Ownership stays with the semantic packages; the façade only
  re-exports.
- **3.2 `multivar.advanced`.** *Done* for the opaque spectral results:
  `typedFrame` and `svdResult` on `PcaFit` and `SvdFit`. Further escape hatches
  land with the remaining Phase 2 opaque results.
- **3.3 `multivar.syntax.unsafe`.** *Done.* `orThrow` and `MultivarException`.
- **3.4 Hierarchy policy.** *Done.* `multivar-package-hierarchy.md` records the
  curated façade, forbids a wildcard root mirror, and adds `analysis`,
  `advanced`, and `syntax.unsafe` to the ownership table.
- **3.5 Ordinary documentation.** *Done.* Ordinary `site-docs/` examples and
  the README import `multivar.analysis.*`. The typed/advanced guide keeps
  explicit `ComponentCount` / `MatrixView` imports beside the façade.

Phase 2 items 2.1–2.7 are landed. The façade re-exports the opaque ordinary
fits; advanced escape hatches open typed frames and operator payloads without
putting them on the ordinary surface.

### Phase 4 — Release topology — **landed for local/CI**

- **4.1 Resolve the Gale dependency.** *Done for the build graph.* Multivar
  depends on `io.github.canardlapin:::gale-core` at a revision-qualified
  version (`1.0.0-<sha12>`) so published POMs declare a Maven coordinate
  rather than a Git `ProjectRef`. Until Gale is on Maven Central,
  `tools/publish-gale-local.sh` installs that exact revision locally. CI runs
  the script before every compile. Switching the same coordinate to Central is
  a Gale release step, not another multivar graph change.
- **4.2 Consumer smoke project.** *Done.* `modules/smoke` depends only on
  publishedLocal `multivar-core` and `multivar-ir` (no in-repo `dependsOn`).
  `sbt smokeCheck` publishLocals the JVM artifacts and compiles documentation
  examples against that graph; CI runs it after the main test matrix.
- **4.3 Binary compatibility.** *Scaffolded.* sbt-mima is wired with empty
  `mimaPreviousArtifacts` and `mimaFailOnNoPrevious := false` until `0.1.0`
  exists. `sbt mimaCheck` is a CI no-op for now; the Phase 0 public-surface
  snapshot remains the pre-release gate. After the first release, point
  `mimaPreviousArtifacts` at `0.1.0`.

Size: medium. Central publication of Gale remains outside this repository and
is the remaining blocker for a consumer who has never run
`publish-gale-local.sh`.

### Phase 5 — PLS regression as the architectural test

Landed. Locality held: `PreparedPair` is shared by PLSC, CCA, RRR, and PLS;
SIMPLS lives in `family.paired` without a new `PairedProgramKind`; façade and
workflow each gained one registration; R parity and law tests live beside the
existing paired suites.

- **5.1 Extract `PreparedPair`** as a family-private value holding the original
  and working views, both preprocessors, the row relationship, moments, and
  schemas. Branch it into the existing cross-spectral problem (PLSC, CCA, RRR)
  and a new predictive latent problem. Preprocessing, row alignment, schemas,
  cross-products, provenance, and raw-coordinate restoration are shared; solver
  and component semantics are not.
- **5.2 Implement SIMPLS** for complete dense multivariate responses, with
  NIPALS deferred so an iterative loop is not part of the first public
  contract. PLS regression must not be a new `PairedProgramKind` case: the
  current PLSC is a one-shot cross-covariance decomposition, structurally
  closer to `PLSSVD` than to PLS regression, whose definition involves
  successive components and response-oriented deflation.
- **5.3 Laws, not examples.** `transform(trainingX)` equals stored `xScores`;
  `predict(x)` equals `x * coefficients + intercept` in raw coordinates; with
  full rank and enough components predictions converge to full OLS; component
  ordering is deterministic; sign orientation is canonical or documented as
  indeterminate; score orthogonality and deflation residuals satisfy the
  algorithm's equations; single- and multi-response cases agree where they
  should; an independent implementation matches predictions, scores up to sign,
  and coefficients. The existing R parity fixtures under `tools/r-parity/` are
  the model to follow.
- **5.4 Measure the non-locality.** A successful addition should require local
  changes to `family.paired`, one façade export, one workflow registration, and
  tests. If it instead requires edits scattered across unrelated layers, that
  is the signal to do Phase 6 before adding more algorithms.

### Phase 6 — Open runtime extension versus closed portable IR

Deferred, and smaller than the review implies. The review counts the IR among
the layers an estimator must touch; it does not, because operator plans
deliberately do not serialize the estimator enum. The real cost is
concentrated in `workflow/Plans.scala`, which has roughly eight to ten
exhaustive `match` sites over `MultivarEstimator` and `FitArtifactKind`.

Separate the two concepts: an open runtime estimator trait with an associated
`Fit` type, so a custom estimator can run locally without automatically having
an IR representation; and a closed, versioned, serializable built-in request
enum for portable workflows, entered by supplying an explicit codec and
executor registration. `MathematicalModelFamily` stays closed and useful as a
curated theorem and guarantee taxonomy, and stops being the gate every runtime
extension passes through.

### Phase 7 — Task façades for advanced families

Deferred. Structured factors currently require an `OpTable`, two geometries,
two penalty structures, and an initialization before solving; GLRM requires
assembling every semantic component by hand. Both should acquire task builders
that infer routine identities, initialization, default loss scaling, and
storage choices while keeping estimand-changing decisions explicit. The
explicit programs remain underneath. This is the intended growth pattern for
the project: establish the exact algebra first, then add a façade that hides
identities and boilerplate but never an estimand.

## Release gate for 0.1

Phases 0 through 4 are the gate. Phase 5 onward is post-0.1.

1. No ordinary fitted result has a public constructor, `copy`, or `unapply`.
2. No ordinary signature mentions `SemanticSpace`, `Op`, `OperatorProgramFit`,
   `SvdResult`, certificates, or path-dependent core types, enforced by the
   Phase 0 surface guard.
3. RRR satisfies the raw-coordinate prediction law under every invertible
   affine preprocessing scheme, and the documentation matches.
4. GPCA centering is explicit, and a nonuniform row metric cannot silently
   receive ordinary centering.
5. Nyström neither discards nor reorders a landmark request, and kernel
   asymmetry beyond tolerance is rejected rather than averaged away.
6. Response-side preprocessing is invertible at fit time.
7. Vocabulary is normalized, and `multivar.analysis` plus `multivar.advanced`
   exist with the documentation examples importing only the façade.
8. `sbt compileAll testAll` and `sbt docsCheck` pass on JVM and Scala.js, the
   signature snapshot is current, and the consumer smoke project compiles
   against published artifacts.

## Decisions taken

1. **`multivar.analysis` is adopted** as the curated façade, together with the
   amendment to the import policy in `multivar-package-hierarchy.md` that
   permits a named-selector façade while continuing to forbid a wildcard root
   mirror.
2. **GPCA `Auto` centering refuses to guess.** It resolves to ordinary
   centering under an identity row metric and returns a typed error under a
   nonuniform row metric, rather than dispatching on metric shape. A caller
   with a nonuniform row metric must name the centering they want.
3. **`Lda` is retired outright.** `FisherDiscriminant` owns the projection and
   `LdaClassifier` owns the classifier. The name `Lda` is never reused, so no
   existing code or documentation can silently change meaning.
4. **CPCA keeps one name.** `Cpca.fit` gains a `preprocessing` parameter
   defaulting to `Center`; the as-supplied behaviour remains available as
   `preprocessing = PreprocessSpec.Pass`. No separate `ConstrainedSvd` estimand
   is introduced.

## Open questions

1. **0.1 scope.** Whether `LdaClassifier`, the estimator-instance form
   (`Pca(components = 5).fit(x)`), and PLS regression are inside the 0.1 gate
   or follow it. All three are additive rather than breaking, so including them
   delays the release without increasing later migration cost; excluding them
   means 0.1 ships without a classifier or a predictive PLS. This document
   currently sequences all three after 0.1 and needs amending if they are
   pulled forward.
2. **Gale release topology**, the one item with a dependency outside this
   repository. It is sequenced late but should be started early, because
   multivar cannot publish a resolvable POM until Gale's own release path
   exists.
