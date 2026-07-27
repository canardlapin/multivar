# Errors and validation

Public smart constructors and fitting methods report invalid requests with a
typed `Left`. Numerical kernels are called only after the public boundary has
checked the conditions it can know.

## Common checks

| Condition | Typical error |
| --- | --- |
| zero or negative dimension | `InvalidDimension` |
| zero components or components above a rank bound | `InvalidComponentRequest` |
| unequal paired row counts | `MatrixShapeMismatch` |
| row or feature metric with the wrong dimension | `MetricShapeMismatch` |
| negative ridge or spectral tolerance | `InvalidTolerance` |
| duplicate or out-of-range landmark | landmark or index error |
| new data with the wrong feature count | `MatrixShapeMismatch` |
| forbidden dense materialization | `DensificationRejected` |
| non-finite matrix entries | finite-value validation error |
| failed eigensystem or decomposition | `SolverFailed` |

## Handle the value at the useful level

Use `map` when only a successful field is needed:

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.analysis.*

val x = Matrix(3, 2)(
  1.0, 0.0,
  0.0, 1.0,
  1.0, 1.0
)

val scoreCount = Pca.fit(x, 1).map(_.scores.cols)
```

```scala mdoc
scoreCount
```

Use `for` when the next operation can fail:

```scala mdoc:silent
val projected =
  for
    fit <- Pca.fit(x, 1)
    scores <- fit.transform(x)
  yield scores
```

```scala mdoc
projected.map(scores => (scores.rows, scores.cols))
```

Match on the result when the program must choose a recovery action:

```scala mdoc
Pca.fit(x, 0) match
  case Right(_)    => "fitted"
  case Left(error) => error.message
```

An iterative fit may be numerically valid without having converged to the
requested tolerance. Inspect its stopping reason and achieved guarantee rather
than treating every `Right` as a global optimum.
