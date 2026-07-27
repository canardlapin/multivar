# Generalized and constrained PCA

Ordinary PCA uses Euclidean geometry. GPCA changes the row and feature
geometries. CPCA restricts the result to supplied row and feature design
spaces.

## Generalized PCA

For diagonal geometries, pass one positive weight per row and per feature:

```scala mdoc:silent
import gale.linalg.{DVec, Matrix}
import multivar.analysis.*

val x = Matrix(4, 3)(
   1.0,  2.0, 0.0,
   3.0, -1.0, 1.0,
   0.5,  4.0, 2.0,
  -2.0,  1.5, 3.0
)

val rowWeights = DVec.fromSeq(Seq(2.0, 1.0, 3.0, 0.5))
val featureWeights = DVec.fromSeq(Seq(1.0, 2.0, 0.75))

val gpca =
  Gpca.fit(
    input = x,
    components = 2,
    rowWeights = rowWeights,
    featureWeights = featureWeights,
    centering = GpcaCentering.ByRowMeasure
  )
```

```scala mdoc
gpca.map(fit => (fit.scores.cols, fit.weights.cols, fit.eigenvalues.length))
```

The dense entry point also accepts `MetricSpec` values for diagonal, dense
symmetric, or operator-backed metrics. A row metric must match the number of
rows; a feature metric must match the number of columns.

`GpcaCentering.Auto` (the default) applies ordinary centering only under an
identity row metric and returns a typed error otherwise, so a nonuniform row
geometry cannot silently change the estimand. Choose `Ordinary`,
`ByRowMeasure`, `OrthogonalToConstant`, `None`, or `AlreadyCentered`
explicitly when the row metric is not identity.

Use `SemanticGpca` when named row and feature spaces, centering evidence, or a
policy for singular geometry is part of the model. Dense GPCA supplies local
identities and delegates to the same checked problem.

## Constrained PCA

CPCA takes a row design and a feature design. Their column spaces define the
allowed row and feature directions.

```scala mdoc:silent
val rowDesign = Matrix(4, 2)(
  1.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
  0.0, 0.0
)

val featureDesign = Matrix(3, 2)(
  1.0, 0.0,
  0.0, 0.0,
  0.0, 1.0
)

val cpca =
  Cpca.fit(
    input = x,
    rowDesign = rowDesign,
    featureDesign = featureDesign,
    components = 2
  )
```

```scala mdoc
cpca.map(fit => (fit.scores.cols, fit.loadings.cols, fit.singularValues.length))
```

`Cpca.fit` centres by default (`PreprocessSpec.Center`). It returns the
complete constrained block. Call `fitBlocks` when the projected and residual
row-by-feature blocks must be estimated separately. Reconstruction vocabulary
distinguishes working-space metric coordinates (`reconstructWorking`) from
original feature coordinates after inverse preprocessing (`reconstruct`), and
whitened versus metric coordinates on each block result.
