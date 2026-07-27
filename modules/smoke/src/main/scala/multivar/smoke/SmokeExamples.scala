package multivar
package smoke

import gale.linalg.Matrix
import multivar.analysis.*
import multivar.ir.SchemaVersion

/** Compile-only consumer of published multivar artifacts.
  *
  * This module must not be wired with `dependsOn` against the in-repo source
  * projects. It exists to prove that Ivy/Maven-published `multivar-core` and
  * `multivar-ir` resolve, including their Gale transitive dependency, and that
  * the ordinary documentation examples still typecheck against that graph.
  */
object SmokeExamples:

  def pcaExample: Either[MultivarError, (Int, Int)] =
    val data = Matrix(4, 3)(
      2.5, 2.4, 0.5,
      0.5, 0.7, -0.1,
      2.2, 2.9, 0.8,
      1.9, 2.2, 0.3
    )
    Pca.fit(data, components = 2).map(fit => (fit.scores.rows, fit.loadings.cols))

  def fisherExample: Either[MultivarError, Int] =
    val x = Matrix(6, 2)(
      0.0, 0.0,
      0.2, 0.1,
      -0.1, 0.2,
      2.0, 2.0,
      2.1, 1.9,
      1.9, 2.2
    )
    FisherDiscriminant.fit(x, Vector(0, 0, 0, 1, 1, 1), components = 1).map(_.criterionValues.length)

  def irSchemaSurvivesPublication: SchemaVersion =
    SchemaVersion.v0_1
