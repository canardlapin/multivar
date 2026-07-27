# Nyström kernel approximation

Nyström approximation replaces a full kernel eigensystem with one built from
selected landmark rows. Use it when a kernel representation is useful but the
full square kernel matrix is too large or unnecessary.

## Fit a linear approximation

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.analysis.*

val x = Matrix(6, 2)(
  0.0, 0.0,
  1.0, 0.0,
  0.0, 1.0,
  1.0, 1.0,
  2.0, 1.0,
  1.0, 2.0
)

val approximation =
  Nystrom.fit(
    input = x,
    components = 2,
    landmarks = Vector(0, 2, 4)
  )
```

```scala mdoc
approximation.map(fit => (fit.scores.rows, fit.scores.cols, fit.eigenvalues.length))
```

The landmark indices must be distinct and within the input row range. The
component count cannot exceed the number of landmarks.

`transform(newX)` computes kernels between new rows and the stored landmarks,
then applies the fitted extension weights:

```scala mdoc:silent
val transformed =
  approximation.flatMap(_.transform(Matrix(1, 2)(0.5, 0.5)))
```

```scala mdoc
transformed.map(values => (values.rows, values.cols))
```

## Use an RBF kernel

`fitRbf` checks `gamma` before constructing the kernel:

```scala mdoc:silent
val rbf =
  Nystrom.fitRbf(
    input = x,
    components = 2,
    landmarks = Vector(0, 2, 4),
    gamma = 0.2
  )
```

```scala mdoc
rbf.isRight
```

The current implementation can preprocess input columns, but it does not
double-center a nonlinear kernel matrix. Centering the inputs is not the same
operation as centering an RBF kernel. The fit records this distinction in its
centering diagnostics.
