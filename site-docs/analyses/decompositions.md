# SVD and PCA

Use SVD to decompose a matrix as supplied. Use PCA when the columns should be
centered before the decomposition.

| Method | Default preprocessing | Typical question |
| --- | --- | --- |
| `Svd.fit` | none | What low-rank directions summarize this matrix? |
| `Pca.fit` | column centering | What directions explain variation around the column means? |

## Fit both methods

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.family.spectral.{Pca, Svd}

val x = Matrix(5, 3)(
  1.0, 2.0, 0.0,
  2.0, 1.0, 1.0,
  3.0, 4.0, 1.0,
  4.0, 3.0, 2.0,
  5.0, 5.0, 3.0
)

val svd = Svd.fit(x, components = 2)
val pca = Pca.fit(x, components = 2)
```

```scala mdoc
(svd.map(_.scores.cols), pca.map(_.scores.cols))
```

Both fits expose:

- `scores`: the training rows in component coordinates;
- `loadings`: one component direction per column; and
- `singularValues`: the retained singular values.

Call `project(newX)` to score new rows. The method checks the number of columns
and applies the preprocessor stored during fitting.

## Choose preprocessing deliberately

Pass `PreprocessSpec.Standardize` when variables measured on different scales
should each have unit sample variance. Pass `PreprocessSpec.Pass` to PCA only
when the matrix has already been prepared and that fact is part of the calling
code.

The fitted loadings are directions, not uniquely signed vectors. A solver may
return the negative of a loading and the negative of its scores without
changing the fitted subspace.

Next, read [PLSC, CCA, and reduced-rank regression](paired.md) for analyses
with two matrices measured on the same rows.
