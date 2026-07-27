# Generalized and constrained PCA

Ordinary PCA uses Euclidean geometry. GPCA changes the row and feature
geometries. CPCA restricts the result to supplied row and feature design
spaces.

## Generalized PCA

For diagonal geometries, pass one positive weight per row and per feature:

```scala mdoc:silent
import gale.linalg.{DVec, Matrix}
import multivar.core.PreprocessSpec
import multivar.family.spectral.Gpca

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
    featureWeights = featureWeights
  )
```

```scala mdoc
gpca.map(fit => (fit.scores.cols, fit.weights.cols, fit.eigenvalues.length))
```

The dense entry point also accepts `MetricSpec` values for diagonal, dense
symmetric, or operator-backed metrics. A row metric must match the number of
rows; a feature metric must match the number of columns.

Use `SemanticGpca` when named row and feature spaces, centering evidence, or a
policy for singular geometry is part of the model. Dense GPCA supplies local
identities and delegates to the same checked problem.

## Constrained PCA

CPCA takes a row design and a feature design. Their column spaces define the
allowed row and feature directions.

```scala mdoc:silent
import multivar.family.cpca.Cpca

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

`Cpca.fit` returns the complete constrained block. Call `fitBlocks` when the
projected and residual row-by-feature blocks must be estimated separately.
Each requested block records its inertia and can be reconstructed in the
original coordinates.
