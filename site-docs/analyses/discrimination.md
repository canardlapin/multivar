# Discriminant and canonical methods

Use `Lda.fit` when each row has a class label and the goal is a low-dimensional
linear separation of those classes.

## Fit linear discriminant analysis

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.family.canonical.Lda

val x = Matrix(6, 2)(
   0.0, 0.0,
   0.2, 0.1,
  -0.1, 0.2,
   2.0, 2.0,
   2.1, 1.9,
   1.9, 2.2
)

val labels = Vector(0, 0, 0, 1, 1, 1)
val lda = Lda.fit(x, labels, components = 1)
```

```scala mdoc
lda.map(fit => (fit.scores.rows, fit.weights.rows, fit.criterionValues.length))
```

`scores` contains the training rows in discriminant coordinates. `weights`
contains the feature directions, and `criterionValues` reports the fitted
discriminant roots. `project(newX)` scores new rows after checking their column
count.

The default trace-scaled ridge keeps the within-class scatter invertible in
many small-sample or high-dimensional problems. Pass `ridge = 0.0` to require
an unregularized positive-definite within-class scatter. The fit then returns a
typed error if that condition does not hold.

The maximum useful component count is bounded by both the number of features
and the number of classes minus one.

## When the input is an effect operator

The canonical package also contains three operator-level model families:

- `CanonicalEffectProblem` separates a declared effect operator from a
  residual operator and fits canonical effect directions.
- `CanonicalRootSpectrum` and `ManovaStatistics` summarize fitted canonical
  roots with Wilks, Pillai, Hotelling--Lawley, or Roy statistics.
- `ConstrainedCanonicalProblem` adds explicit frame constraints to a canonical
  effect problem.

These APIs start after a caller has defined the effect, residual, and any
constraints. They do not infer a design matrix or statistical contrast from
labels. Use them when those operators already come from a model-fitting layer;
use `Lda.fit` for labelled dense data.
