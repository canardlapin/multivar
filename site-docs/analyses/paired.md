# PLSC, CCA, and reduced-rank regression

These methods take two matrices with the same number of rows. They differ in
the quantity they optimize and in whether the result is symmetric or
predictive.

| Method | Use it to | Main fitted quantity |
| --- | --- | --- |
| PLSC | find paired directions with high cross-covariance | paired weights and covariances |
| CCA | find paired directions with high correlation | paired weights and correlations |
| Reduced-rank regression | predict `y` from `x` with a rank limit | coefficient matrix |

## Fit PLSC and CCA

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.family.paired.{Cca, Plsc, ReducedRankRegression}

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

Both fits expose scores and weights for each input. `projectX` and `projectY`
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

The fit exposes `coefficients` in the original predictor and response
coordinates. `predict(newX)` applies the fitted predictor preprocessing,
coefficient map, and inverse response preprocessing. Pass `ridge = 1e-3`, for
example, when ordinary least squares is unstable.

All three methods reject unequal row counts. CCA rejects negative ridge values.
Reduced-rank regression also rejects a requested rank above the predictor or
response rank bound.
