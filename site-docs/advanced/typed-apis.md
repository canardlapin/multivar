# Typed and operator-backed APIs

The dense overloads are the normal entry point. Use the typed overloads when a
request will be stored, composed, or evaluated over sparse or lazy data.

## Opt into checked request values

```scala mdoc:silent
import gale.linalg.Matrix
import multivar.analysis.*
import multivar.core.{ComponentCount, MatrixView}

val x = Matrix(4, 2)(
  1.0, 2.0,
  2.0, 1.0,
  3.0, 4.0,
  4.0, 3.0
)

val fit =
  for
    components <- ComponentCount(1)
    pca <- Pca.fit(MatrixView.dense(x), components)
  yield pca
```

```scala mdoc
fit.map(_.scores.cols)
```

The result is the same `PcaFit` returned by the dense overload. The difference
is when validation occurs: a `ComponentCount` can be constructed once and
reused in a workflow.

## Use `MatrixView` for storage-aware work

`MatrixView` represents dense, sparse, selected, transposed, and lazily
preprocessed matrices behind one checked interface. Operations such as
right-multiplication can remain representation-aware. An operation that must
produce a dense matrix accepts a `StoragePolicy`; `RejectDense` turns that
materialization into a typed error.

Named semantic tables add row and feature spaces to the same numerical data.
Use them when feature identity or operator orientation must survive
serialization and composition. Ordinary dense users do not need to create
those spaces.

## Use model specifications for reusable workflows

`ModelSpec` and the plan types in `multivar.workflow` describe a request before
data are executed. They are useful for fold fitting, remote execution, or IR
serialization. Direct `fit` calls are clearer when the data are already in
memory and no reusable plan is needed.
