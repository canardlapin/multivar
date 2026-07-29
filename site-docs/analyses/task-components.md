# Task components

Use `TaskComponents` when a least-squares prediction or table-decomposition
task should be resolved into a scientifically specified focus and its
orthogonal remainder, with exact in-sample loss accounting.

If there is no external scientific resolution, prefer
`ReducedRankRegression.fit` for prediction or `Gpca.fit` for metric
decomposition. `TaskComponents` is not a replacement for those methods, and it
is not a second spelling of `Cpca`. Public `TaskComponents.decompose` requires
a nontrivial resolution; unresolved metric SVD stays with `Gpca`.

| Method | Object | External structure |
| --- | --- | --- |
| `Gpca` | a data table | none |
| `Cpca` | a data table | crossed row × feature blocks |
| `ReducedRankRegression` | a directed prediction task | none |
| `TaskComponents.fit` | a prediction task \(X\to Y\) | one scientific focus plus its orthogonal remainder |
| `TaskComponents.decompose` | a metric table task | one scientific focus plus its orthogonal remainder |

`Cpca` partitions table inertia into crossed `G×H` blocks. `TaskComponents`
resolves a least-squares task operator according to weight, structure,
covariance, or regression semantics. Neither replaces the other.

## Fit a resolved prediction

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.analysis.*

val x = Matrix(6, 2)(
  1.0, 0.0,
  2.0, 1.0,
  3.0, 1.0,
  4.0, 2.0,
  5.0, 3.0,
  6.0, 4.0
)

val y = Matrix(6, 2)(
  2.0, 1.0,
  4.0, 1.5,
  6.0, 2.0,
  8.0, 3.0,
 10.0, 4.0,
 12.0, 5.0
)

val responseGroups = Matrix(2, 1)(
  1.0,
  0.0
)

val fitted =
  TaskComponents.fit(
    predictors = x,
    responses = y,
    components = 1,
    resolution = TaskResolution.Target.weights(responseGroups)
  )
```

```scala mdoc
fitted.map(fit => (fit.focus.rank, fit.remainder.rank, fit.value.availableGain > 0.0))
```

`components = k` means `RankBudget.Focus(k)`: retain focus components only.
Use `RankBudget.Total(k)` or `RankBudget.Split(focus, remainder)` when the
rank plan must be stated explicitly.

A structure resolution differs by one word and optionally names a dual policy:

```scala mdoc:silent
val structured =
  TaskComponents.fit(
    predictors = x,
    responses = y,
    components = 1,
    resolution = TaskResolution.Target.structure(
      responseGroups,
      dual = DualSpectrum.Ridge(1e-4)
    )
  )
```

```scala mdoc
structured.map(_.dualSpectrum.map(_.label))
```

## Source-side resolution

Axis meaning depends on the task:

| Task | Source axis | Target axis |
| --- | --- | --- |
| Predict \(X\to Y\) | predictor features (\(p\)) | response features (\(q\)) |
| Decompose \(X\) | observations (\(n\)) | features (\(p\)) |

For prediction, a source design therefore has one row per predictor feature,
not per observation.

```scala mdoc:silent
val predictorGroups = Matrix(2, 1)(
  1.0,
  0.0
)

val sourceResolved =
  TaskComponents.fit(
    predictors = x,
    responses = y,
    components = 1,
    resolution = TaskResolution.Source.regression(predictorGroups)
  )
```

```scala mdoc
sourceResolved.map(fit => (fit.sourceAxis, fit.taskSupportRank, fit.focus.value.availableGain >= 0.0))
```

`Source.covariance` weights source information by task strength.
`Source.regression` retains only the part of the source design that lies in the
numerical support of \(K\). Passing an \(n\times d\) observation design to a
prediction fit is a typed shape error that names the predictor-feature axis.

## Decompose a resolved table

```scala mdoc:silent
val table = Matrix(6, 3)(
  1.0, 0.0, 2.0,
  2.0, 1.0, 1.0,
  3.0, 1.0, 0.0,
  4.0, 2.0, 1.0,
  5.0, 3.0, 2.0,
  6.0, 4.0, 3.0
)

val featureGroups = Matrix(3, 1)(
  1.0,
  0.0,
  0.0
)

val decomposed =
  TaskComponents.decompose(
    input = table,
    resolution = TaskResolution.Target.weights(featureGroups),
    components = 1
  )
```

```scala mdoc
decomposed.map { fit =>
  (
    fit.sourceAxis,
    fit.targetAxis,
    fit.resolutionTransfer.label,
    fit.transform(table).map(_.cols)
  )
}
```

New observations use the fitted feature-side transfer \(F = Z^{+}EZ\), folded
into original-feature analysis/synthesis maps:

```text
transform(newX)                 // retained component scores
reconstruct(newX)               // baseline + focus + remainder
focusContribution(newX)
remainderContribution(newX)
residual(newX)                  // newX - reconstruct(newX)
focus.transform(newX)           // branch scores
validationLoss(newX)            // held-out Frobenius evaluation (not TaskValue)
resolve(newX)                   // scientific focus / remainder / unsupported
unsupportedFraction(newX)       // whitened energy outside training task support
```

```scala mdoc:silent
val heldOut = Matrix(2, 3)(
  2.5, 1.5, 1.0,
  3.5, 2.0, 2.5
)
```

```scala mdoc
decomposed.flatMap { fit =>
  for
    hat <- fit.reconstruct(heldOut)
    focusScores <- fit.focus.transform(heldOut)
    report <- fit.validationLoss(heldOut)
    parts <- fit.resolve(heldOut)
    frac <- fit.unsupportedFraction(heldOut)
  yield
    (
      hat.rows,
      focusScores.cols,
      report.fittedLoss <= report.baselineLoss,
      parts.unsupported.cols,
      frac >= 0.0
    )
}
```

`resolve` uses the full scientific transfers \(F_f,F_r\) and
\(F_0 = I - P_{Z^\top}\). It is a diagnostic split, not a third estimand branch:
rank-limited `reconstruct` may discard focus/remainder components that `resolve`
still attributes to the scientific parts.

Reconstruction cancels the fitted centering baseline:

```text
reconstruct(x) = baseline(rows) + focusContribution(x) + remainderContribution(x)
```

For decomposition, `Source.*` designs have one row per observation. New-data
projection still applies; `resolutionTransfer` is
`TrainingSourceInduced` and means the feature map learned from the training
source design, not re-evaluation of a new source-design row. Unresolved metric
decomposition remains `Gpca.fit`; `TaskResolution.Whole` is rejected on the
ordinary `decompose` entry point.

Training quantities such as `availableGain` stay in-sample and are not claimed
for held-out rows. Use `validationLoss` for declared out-of-sample reconstruction
evaluation.

## Inspect task gain and predictions

`availableGain` is the in-sample reduction in the declared quadratic loss that
the unrestricted task can achieve. `retainedGain` is the reduction achieved by
the rank-limited fit. These quantities are training loss reductions, not
held-out predictive utility.

```scala mdoc
fitted.map { fit =>
  (
    fit.value.baselineLoss > fit.value.fittedLoss,
    fit.focus.value.shareOfAvailable >= 0.0,
    fit.predict(x).map(values => (values.rows, values.cols))
  )
}
```

Branch contributions cancel the response affine baseline:

```text
predict(x) = baseline + focusContribution(x) + remainderContribution(x)
```

## When to stay with ReducedRankRegression

Under `TaskResolution.Whole`, `TaskComponents.fit` recovers the same predictive
map as `ReducedRankRegression` for matching preprocessing and ridge. The
ordinary chooser still points unresolved prediction to `ReducedRankRegression`
so the direct methods remain the short path.
