# Perturbation inference

The `multivar-inference` artifact adds typed resampling designs, ordered-root
tests, stability summaries, and deterministic Monte Carlo receipts. It is
separate from `multivar-core`: fitting a PCA or CCA model does not implicitly
claim inferential validity.

Add the artifact for your platform:

```scala
libraryDependencies +=
  "io.github.canardlapin" %% "multivar-inference" % multivarVersion
```

A plan names the fitted family, target, null hypothesis, sampling design,
latent-unit policy, Monte Carlo policy, requested evidence, and root seed:

```scala mdoc
import multivar.inference.*
import resample4s.kernel.Seed

val rows = RowCount(24).toOption.get
val draws = MonteCarloDraws(999).toOption.get
val alpha = Alpha(0.05).toOption.get
val batch = BatchSize(25).toOption.get

val spec = InferenceSpec.of(
  FitDescriptor.Pca,
  TargetSpec.VarianceRoots,
  NullSpec.PermuteRows,
  ResamplingDesign.exchangeableRows(rows),
  UnitPolicy.GroupNearTies(RelativeGap(0.01).toOption.get),
  MonteCarloPolicy.Sequential(
    draws,
    alpha,
    batch,
    SequentialBoundary.FromAlpha
  ),
  RequestedEvidence.SignificanceOnly,
  Seed.fromLong(1729L)
)

val program = InferenceCompiler.compile(spec)
program.summary
```

`InferenceSpec` is plan data: it contains no matrices, fitted objects,
callbacks, open resources, or scheduler handles. Compilation rejects
unsupported target, null, and design relationships through typed evidence.

`InferenceExecutor.runProgram` is the significance-ladder interpreter. It
accepts `SignificanceOnly` plans and returns a typed `UnsupportedEvidence`
refusal for stability-only or combined requests. Stability has its own
bootstrap and alignment operations; compiling a request does not imply that
every interpreter can satisfy it.

The current significance protocols define axis statistics. A
`GroupNearTies` or declared policy remains valid plan data, but if it actually
forms a multi-axis subspace the executor returns `UnsupportedUnitPolicy`.
Subspace significance needs a rotation-invariant statistic; Multivar does not
attach axis-wise p-values to an unoriented subspace.

Execution uses Resample4s seed paths, so each action is determined by the root
seed and `ReplicateId`, independent of task order. Receipts record both the
randomization algorithm and seed-derivation algorithm.

Conditioned designs name a projector/whitening resource but do not smuggle that
runtime object into the plan. Prepare those resources explicitly with
`ResidualPermutationAction`; the ordinary ladder executor refuses a
conditioned program if they have not been supplied.

Bootstrap stability remains separate from significance testing. Selection
frequency, aligned loadings, and principal angles are stability evidence and
are not converted into ordinary p-values.
