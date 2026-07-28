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
import multivar.analysis.*

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
- `loadings`: one component direction per column;
- `components`: the same directions transposed, one per row;
- `singularValues`: the retained singular values;
- `center` and `scale`: what the fitted preprocessing removed, as
  `(x - center) / scale`; and
- `nFeatures`, `requestedComponents`, `effectiveComponents`.

Call `transform(newX)` to score new rows. The method checks the number of
columns and applies the preprocessor stored during fitting.
`inverseTransform(scores)` maps component coordinates back to original feature
coordinates, and `reconstruct(newX)` composes the two to give the rank-`k`
approximation of `newX`.

```scala mdoc:silent
val approximation =
  pca.flatMap(_.reconstruct(x))
```

```scala mdoc
approximation.map(values => (values.rows, values.cols))
```

The approximation lives in the same units as the input matrix, including when
preprocessing was `Standardize()` rather than centering alone.

The two fits differ in what they claim about the retained shares. PCA centres,
so it reports `explainedVariance` and `explainedVarianceRatio`. SVD does not
centre, so its shares are sums of squares about the origin rather than
variances, and it reports `inertia` and `inertiaRatio` instead. Both describe
the preprocessed data: under `Standardize()` the PCA variances are those of the
standardized columns, which is the correlation-PCA convention. The ratios use
the total of the full preprocessed matrix as their denominator, so a truncated
fit's shares sum to less than one rather than being renormalized.

Because a variance cannot be estimated from a single observation, `Pca.fit`
requires at least two rows. `Svd.fit`, which makes no variance claim, does not.

## Choose preprocessing deliberately

Pass `PreprocessSpec.Standardize()` when variables measured on different scales
should each have unit sample variance. Pass `PreprocessSpec.Pass` to PCA only
when the matrix has already been prepared and that fact is part of the calling
code.

The fitted loadings are directions, not uniquely signed vectors. A solver may
return the negative of a loading and the negative of its scores without
changing the fitted subspace.

Next, read [PLSC, CCA, RRR, and PLS](paired.md) for analyses
with two matrices measured on the same rows.
