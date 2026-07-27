package multivar
package advanced

import gale.linalg.Matrix
import multivar.analysis.*
import multivar.core.{MatrixView, MultivarError}
import multivar.family.canonical.{TraceRidgeFraction, WithinScatterPolicy}
import multivar.family.cpca.CpcaFit
import multivar.syntax.unsafe.orThrow

import gale.linalg.DMat

class AdvancedSurfaceSuite extends munit.FunSuite:

  private val data = Matrix(4, 3)(
    2.5, 2.4, 0.5,
    0.5, 0.7, -0.1,
    2.2, 2.9, 0.8,
    1.9, 2.2, 0.3
  )

  private val pairedY = Matrix(4, 2)(
    1.0, 0.5,
    0.2, -0.3,
    1.1, 0.7,
    0.8, 0.4
  )

  private def assertMatrixClose(actual: DMat, expected: DMat, tol: Double): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assertEqualsDouble(actual(row, col), expected(row, col), tol)
        col += 1
      row += 1

  test("advanced opens the typed frame and SVD payload from opaque spectral fits"):
    val pca = Pca.fit(data, components = 2).orThrow
    val svd = Svd.fit(data, components = 2).orThrow

    assertEquals(pca.typedFrame.diagnostics.method, "pca")
    assertEquals(svd.typedFrame.diagnostics.method, "svd")
    assertEquals(pca.svdResult.singularValues.length, pca.effectiveComponents)
    assertEquals(svd.svdResult.singularValues.length, svd.effectiveComponents)
    assertEquals(pca.typedFrame.scores.rows, data.rows)
    assertEquals(svd.typedFrame.scores.rows, data.rows)

  test("advanced escape hatches agree with ordinary PLS, Fisher, and CPCA surfaces"):
    val pls = PlsRegression.fit(data, pairedY, components = 1).orThrow
    assertMatrixClose(pls.coefficientTransform.predict(data).orThrow, pls.predict(data).orThrow, 1e-12)
    assertEquals(pls.sourceFrame.diagnostics.method, "pls.source")

    val labels = Vector(0, 0, 1, 1)
    val fisher = FisherDiscriminant
      .fit(
        data,
        labels,
        components = 1,
        withinRegularization = WithinScatterPolicy.FixedTraceScaledRidge(TraceRidgeFraction.unsafe(1e-3))
      )
      .orThrow
    assertMatrixClose(fisher.typedFrame.project(data).orThrow, fisher.transform(data).orThrow, 1e-12)

    val rowDesign = Matrix(4, 1)(1.0, 0.0, 0.0, 0.0)
    val featureDesign = Matrix(3, 1)(1.0, 0.0, 0.0)
    val cpca = Cpca.fit(data, rowDesign, featureDesign, components = 1, PreprocessSpec.Center).orThrow
    val restored = cpca.fittedPreprocessor
      .inverseTransform(MatrixView.dense(cpca.reconstructWorking().orThrow))
      .flatMap(_.toDense())
      .orThrow
    assertMatrixClose(cpca.reconstruct().orThrow, restored, 1e-10)
    assertEquals(CpcaFit.preprocessorOf(cpca), cpca.fittedPreprocessor)

  test("orThrow surfaces a MultivarError as MultivarException"):
    val failure =
      try
        Pca.fit(data, components = 0).orThrow
        None
      catch
        case error: multivar.syntax.unsafe.MultivarException => Some(error.error)

    failure match
      case Some(MultivarError.InvalidDimension("component count", 0)) => ()
      case other => fail(s"expected InvalidDimension, got $other")
