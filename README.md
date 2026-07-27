# multivar

Multivar is a Scala 3 library for multivariate analysis on the JVM and
Scala.js. It provides familiar methods such as PCA, PLSC, CCA, and reduced-rank
regression, together with generalized, sparse, kernel, and multiblock models.

[Read the guide](https://canardlapin.github.io/multivar/) or begin with the
short PCA example below.

The library keeps preprocessing, matrix geometry, fitted transformations, and
solver evidence in the result. A fitted model therefore records what was
estimated and checks whether new data have the expected shape. APIs that accept
named feature schemas also check feature identity.

## Status

Multivar is in early development on the 0.1 line. It is not yet published to
Maven Central. The source build uses Scala 3.7.4 and supports the JVM and
Scala.js.

The first release will publish two artifacts. JVM projects will use:

```scala
libraryDependencies += "io.github.canardlapin" %% "multivar-core" % version
libraryDependencies += "io.github.canardlapin" %% "multivar-ir" % version
```

Scala.js cross-projects will use `%%%` in place of `%%`.

## PCA in a few lines

After cloning the repository, start `sbt coreJVM/console` and run:

```scala
import gale.linalg.Matrix
import multivar.family.spectral.Pca

val data = Matrix(4, 3)(
  2.5, 2.4,  0.5,
  0.5, 0.7, -0.1,
  2.2, 2.9,  0.8,
  1.9, 2.2,  0.3
)

val analysis = Pca.fit(data, components = 2)
```

`Pca.fit` centers each column by default. It returns
`Either[MultivarError, PcaFit]`, so an invalid component count or matrix shape
is a value to handle rather than an exception from the numerical kernel. A
successful fit exposes the results directly:

```scala
analysis match
  case Right(pca) =>
    val scores = pca.scores
    val loadings = pca.loadings
    val singularValues = pca.singularValues
  case Left(error) =>
    println(error.message)
```

`scores` contains the training observations in principal-component
coordinates. There is no need to project the training matrix again. To score
new observations, call `pca.project(newData)`; it applies the fitted centering
and loadings.

Use `PreprocessSpec.Pass` to leave the columns unchanged or
`PreprocessSpec.Standardize` to scale them to unit sample variance:

```scala
import multivar.core.PreprocessSpec

val standardized =
  Pca.fit(data, components = 2, preproc = PreprocessSpec.Standardize)
```

## Analyses

| Task | Main entry points |
| --- | --- |
| PCA and SVD | `Pca`, `Svd` |
| Paired latent analysis | `Plsc`, `Cca` |
| Directed prediction | `ReducedRankRegression` |
| Linear discrimination | `Lda` |
| Generalized geometry | `Gpca`, `Cpca` |
| Kernel approximation | `Nystrom` |
| Sparse and smooth factors | `RankOneStructuredFactorization` |
| Mixed data and missing entries | `GeneralizedLowRankProgram` |
| Multiblock analysis | `MultisetAssociation`, `AlignedSharedScoreGlrm` |

See the [guide](site-docs/getting-started.md) for a first analysis and the
[method chooser](site-docs/reference/choosing-method.md) for the rest of the
library.

The direct `DMat` overload is intended for ordinary dense analysis. Sparse and
lazy operator-backed data use the `MatrixView` interface. Operations that may
materialize those inputs expose a `StoragePolicy` at that boundary, so callers
can allow materialization or require the operation to fail.

## Results and errors

The main smart constructors and fitting methods validate dimensions, component
counts, regularization values, and matrix roles. Feature-aware operations also
validate feature identities. These methods return typed errors that describe
the failed condition.

Fitted transformations retain their preprocessing and training schema. SVD,
PCA, PLSC, CCA, LDA, and GPCA fits can score new observations. Nyström fits
transform new observations through their stored landmarks. Reduced-rank
regression fits predict responses. Reconstruction and nonlinear encoding remain
separate fitted capabilities.

Iterative models report the stopping reason and the guarantee achieved.
Reaching an iteration limit does not become a convergence claim, and
a stationary solution does not become a claim of global optimality.

## Modules

`multivar-core` contains the analysis APIs, typed matrix geometry,
statistical solver adapters, and fitted transformations.

`multivar-ir` contains JSON codecs, schemas, and conformance fixtures for
portable analysis programs and evidence. Most Scala applications need only
`multivar-core`.

Multivar depends on [Gale](https://github.com/canardlapin/gale) for portable
matrices, linear operators, and general numerical solvers. Dataset access,
storage, scheduling, and domain-specific adapters remain outside this library.

## Build

Run the complete JVM and Scala.js build with:

```sh
sbt compileAll testAll
```

Compile the shared API reference, execute the guide examples, and render the
local documentation site with:

```sh
sbt docsCheck
```

Numerical changes should include an analytic law, an adversarial case, or an
independent reference result. See [CONTRIBUTING.md](CONTRIBUTING.md) before
changing public APIs or wire formats.

## License

Licensed under the [Apache License 2.0](LICENSE).
