# PLSC, CCA, reduced-rank regression, and PLS

These methods take two matrices with the same number of rows. They differ in
the quantity they optimize and in whether the result is symmetric or
predictive.

| Method | Use it to | Main fitted quantity |
| --- | --- | --- |
| PLSC | find paired directions with high cross-covariance | paired weights and covariances |
| CCA | find paired directions with high correlation | paired weights and correlations |
| Reduced-rank regression | predict `y` from `x` with a rank limit | coefficient matrix |
| PLS regression (SIMPLS) | predict `y` from `x` through successive latent components | coefficients, X-scores, rotations |

## Fit PLSC and CCA

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

val plsc = Plsc.fit(x, y, components = 1)
val cca = Cca.fit(x, y, components = 1)
```

```scala mdoc
(plsc.map(_.covariances.length), cca.map(_.correlations.length))
```

Both fits expose scores and weights for each input. `transformX` and `transformY`
apply the corresponding fitted transform to new rows.

CCA uses a small symmetric ridge by default. Supply one value to change both
sides, or two values when the covariance estimates need different
regularization:

```scala mdoc:silent
val regularized =
  Cca.fit(x, y, components = 1, xRidge = 1e-4, yRidge = 1e-3)
```

```scala mdoc
regularized.isRight
```

## Predict with reduced-rank regression

Reduced-rank regression gives the inputs different roles: `x` contains
predictors and `y` contains responses.

```scala mdoc:silent
val regression =
  ReducedRankRegression.fit(x, y, components = 1)

val fittedValues =
  regression.flatMap(_.predict(x))
```

```scala mdoc
fittedValues.map(values => (values.rows, values.cols))
```

`coefficients` and `intercept` are in original predictor/response units, so
`predict(newX)` matches `newX * coefficients + intercept` after accounting for
row layout. `workingCoefficients` is the same map in preprocessed coordinates.
Both sides default to `PreprocessSpec.Center`.

Because `predict` must undo the response preprocessing, the fit requires a
response spec whose fitted form is invertible. A response weight of zero is
therefore rejected when the model is fitted rather than when the first
prediction is made.

## Predict with PLS regression

PLS regression builds successive latent components that explain covariance
between `x` and `y`. The first public algorithm is SIMPLS (`PlsAlgorithm.Simpls`);
NIPALS is deferred.

```scala mdoc:silent
val pls =
  PlsRegression.fit(x, y, components = 1)

val plsPredictions =
  pls.flatMap(_.predict(x))

val plsScores =
  pls.flatMap(_.transform(x))
```

```scala mdoc
(
  plsPredictions.map(values => (values.rows, values.cols)),
  plsScores.map(values => (values.rows, values.cols))
)
```

`transform` returns X-scores for new predictors. `coefficients` / `intercept`
are again in original units. Component signs are oriented so the
largest-magnitude entry of each X-rotation is non-negative. The retained
component count is capped by `min(n, p)` (predictor space), so univariate
multi-component fits are allowed.

SIMPLS is checked against R package `pls` via `plsr(..., method = "simpls")`
fixtures under `tools/r-parity/`.

PLS regression is not a one-shot cross-covariance SVD: it is a separate
predictive estimator that shares preparation with PLSC, CCA, and RRR through
`PreparedPair`, but does not add a `PairedProgramKind` case.

## Shared rejection rules

All four methods reject unequal row counts. CCA rejects negative ridge values.
Reduced-rank regression and PLS also reject a requested rank above the
predictor/response rank bound.
