# Getting started

This page fits PCA to a dense matrix, reads the fitted scores, and applies the
fit to new rows.

## Run Multivar from source

Multivar does not yet have a Maven Central release. Clone the repository and
start a JVM console:

```sh
git clone https://github.com/canardlapin/multivar.git
cd multivar
sbt coreJVM/console
```

Applications will use these coordinates after the first release:

```scala
libraryDependencies += "io.github.canardlapin" %% "multivar-core" % version
```

Scala.js cross-projects use `%%%` instead of `%%`.

## Build a data matrix

Rows are observations. Columns are variables.

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.family.spectral.Pca

val data = Matrix(6, 3)(
  2.5, 2.4,  0.5,
  0.5, 0.7, -0.1,
  2.2, 2.9,  0.8,
  1.9, 2.2,  0.3,
  3.1, 3.0,  1.0,
  1.2, 1.1,  0.0
)
```

## Fit PCA

```scala mdoc:silent
val result = Pca.fit(data, components = 2)
```

`Pca.fit` centers each column. It returns
`Either[MultivarError, PcaFit]`. The component count is checked before the
numerical solver runs.

Read a field with `map` when no later operation can fail:

```scala mdoc
val scoreShape = result.map(pca => (pca.scores.rows, pca.scores.cols))
```

Use a `for` expression when the next operation also returns `Either`. Projection
checks that the new matrix has the same number of columns as the training
matrix, then applies the fitted centering and loadings.

```scala mdoc:silent
val newData = Matrix(2, 3)(
  2.0, 2.1, 0.4,
  0.8, 0.9, 0.0
)

val newScores =
  for
    pca <- result
    scores <- pca.project(newData)
  yield scores
```

```scala mdoc
newScores.map(scores => (scores.rows, scores.cols))
```

## Handle an invalid request

Invalid dimensions and component counts are values, not exceptions from a
matrix kernel.

```scala mdoc
Pca.fit(data, components = 0).left.map(_.message)
```

For a complete account of result handling and preprocessing, read
[Data, preprocessing, and fitted results](concepts/data-and-results.md). To
choose another analysis, see the [method chooser](reference/choosing-method.md).
