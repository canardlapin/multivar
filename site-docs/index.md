# Multivar

Multivar is a Scala 3 library for multivariate analysis on the JVM and
Scala.js. It includes PCA, PLSC, CCA, reduced-rank regression, PLS regression,
LDA, generalized and constrained PCA, kernel approximations, structured factors,
generalized low-rank models, and multiblock analysis.

The ordinary dense APIs accept Gale matrices and plain component counts. They
return fitted results whose fields have statistical names such as `scores`,
`loadings`, `correlations`, and `coefficients`.

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.analysis.*

val observations = Matrix(6, 3)(
  2.5, 2.4,  0.5,
  0.5, 0.7, -0.1,
  2.2, 2.9,  0.8,
  1.9, 2.2,  0.3,
  3.1, 3.0,  1.0,
  1.2, 1.1,  0.0
)

val analysis = Pca.fit(observations, components = 2)
```

```scala mdoc
analysis.map(pca => (pca.scores.rows, pca.scores.cols))
```

`Pca.fit` centers the columns and validates the requested rank. On success,
the example has six rows of scores and two principal components. On failure,
it returns a `MultivarError` in `Left`.

Start with [Getting started](getting-started.md). If you already know the
analysis you need, use the [method chooser](reference/choosing-method.md).

## Status

Multivar is early-development research software on the 0.1 line. Its APIs,
numerical behavior, and result types may change without a migration period.
Use it for research and evaluation, and validate results independently before
relying on them in consequential work.

No release has been published to Maven Central. Development builds remain on
the 0.1 line until the API and numerical contracts are stable enough for a
tagged release. The source build uses Scala 3.7.4 and supports the JVM and
Scala.js.

## Main tasks

| Task | Guide |
| --- | --- |
| Reduce one matrix to a few components | [SVD and PCA](analyses/decompositions.md) |
| Relate or predict two matrices | [PLSC, CCA, RRR, and PLS](analyses/paired.md) |
| Separate labelled groups | [Discriminant and canonical methods](analyses/discrimination.md) |
| Supply row or feature geometry | [Generalized and constrained PCA](analyses/geometry.md) |
| Approximate a kernel eigensystem | [Nyström approximation](analyses/kernels.md) |
| Fit sparse or smooth factors | [Structured factors](advanced/structured-factors.md) |
| Model mixed data or missing cells | [Generalized low-rank models](advanced/mixed-data.md) |
| Combine several data blocks | [Multiblock models](advanced/multiblock.md) |
