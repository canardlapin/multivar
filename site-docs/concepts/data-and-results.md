# Data, preprocessing, and fitted results

Most analyses begin with a `gale.linalg.DMat`. Rows are observations and
columns are variables. Dense fitting methods accept this matrix directly.

## Preprocessing belongs to the fit

PCA centers columns by default. SVD and Nyström leave them unchanged. Methods
with two inputs fit a separate preprocessor to each input.

Pass a `PreprocessSpec` when the default is not right for the analysis:

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.analysis.*

val x = Matrix(5, 3)(
  1.0, 10.0, 2.0,
  2.0, 12.0, 1.0,
  3.0, 14.0, 4.0,
  4.0, 16.0, 3.0,
  5.0, 18.0, 5.0
)

val standardized =
  Pca.fit(x, components = 2, preproc = PreprocessSpec.Standardize())
```

`Standardize()` subtracts the training mean and divides by the standard
deviation. Its argument names the denominator: `VarianceConvention.Sample`, the
default, divides the centered sum of squares by `n - 1`, and
`VarianceConvention.Population` divides by `n`. `Center` only subtracts the
mean. `MultiplyColumns` applies weights you supply. `Pass` leaves values
unchanged. The fitted preprocessor is stored with the result, so `transform`,
`project`, and `predict` apply the same operation to new rows.

Applying preprocessing always works; undoing it does not. A zero column weight
transforms fine and cannot be reversed, so the ability to return original
coordinates lives in the type rather than in a runtime check: only a
`FittedInvertiblePreprocessor` has `inverseTransform`, and
`requireInvertible` is the one way to obtain one. Analyses that must report
results in original units, such as the response side of a regression, demand
that type in their own signatures, so an unusable scale is rejected when it is
supplied rather than at the first prediction.

## Results name the quantities users inspect

The direct fits expose their common outputs without requiring a second
projection of the training data:

| Fit | Main fields |
| --- | --- |
| `PcaFit` | `scores`, `loadings`, `singularValues`, `explainedVariance`, `explainedVarianceRatio`, `center`, `scale` |
| `SvdFit` | `scores`, `loadings`, `singularValues`, `inertia`, `inertiaRatio`, `center`, `scale` |
| `PlscFit` | `xScores`, `yScores`, `xWeights`, `yWeights`, `covariances` |
| `CcaFit` | `xScores`, `yScores`, `xWeights`, `yWeights`, `correlations` |
| `ReducedRankRegressionFit` | `coefficients`, `intercept`, scores, weights |
| `PlsRegressionFit` | `coefficients`, `intercept`, X-scores, rotations |
| `FisherDiscriminantFit` | `scores`, `weights`, `criterionValues` |
| `GpcaFit` | `scores`, `weights`, `eigenvalues`, `singularValues`, `center` |
| `CpcaFit` | `scores`, `loadings`, `singularValues` |
| `NystromFit` | `scores`, `eigenvalues` |

The result also keeps diagnostics, provenance, and solver evidence when the
method produces them. Those records distinguish a requested guarantee from
the guarantee the solver achieved.

## Errors compose with the result

Smart constructors and fitting methods return `Either`. Use `map` to read a
successful result, `flatMap` or a `for` expression to continue with another
checked operation, and pattern matching when success and failure need different
actions.

```scala mdoc
standardized match
  case Right(pca) =>
    s"${pca.scores.rows} rows, ${pca.scores.cols} components"
  case Left(error) =>
    error.message
```

The [error reference](../reference/errors.md) lists the checks made at the
public boundary.

## Dense and typed inputs serve different jobs

Use the dense overloads for ordinary in-memory analysis. `MatrixView`,
`ComponentCount`, named spaces, and operator-backed metrics are available when
an analysis must preserve sparse storage, carry feature identity, or be stored
as a reusable program. They are not prerequisites for PCA on a dense matrix.

Read [Typed and operator-backed APIs](../advanced/typed-apis.md) when those
requirements arise.
