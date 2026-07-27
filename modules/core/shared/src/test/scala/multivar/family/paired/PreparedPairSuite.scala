package multivar
package family.paired

import multivar.core.*

import gale.linalg.DMat

/** Contract tests for the shared paired preparation seam. */
class PreparedPairSuite extends munit.FunSuite:

  private val x = GaleNumerics.matrixFromRows(
    Vector(
      Vector(1.0, 2.0),
      Vector(3.0, 4.0),
      Vector(5.0, 6.0),
      Vector(7.0, 8.0)
    )
  )
  private val y = GaleNumerics.matrixFromRows(
    Vector(
      Vector(0.5, -1.0),
      Vector(1.5, 0.0),
      Vector(2.5, 1.0),
      Vector(3.5, 2.0)
    )
  )

  private def view(matrix: DMat): MatrixView =
    MatrixView.dense(matrix)

  test("from and fromPredictive reject unequal rows and metric shape mismatches"):
    val shortY = GaleNumerics.matrixFromRows(Vector(Vector(1.0, 2.0), Vector(3.0, 4.0)))
    PreparedPair.from(view(x), view(shortY), PreprocessSpec.Center, PreprocessSpec.Center) match
      case Left(MultivarError.MatrixShapeMismatch(detail)) =>
        assert(detail.contains("equal rows"), detail)
      case other =>
        fail(s"expected row mismatch, got $other")

    val badMetric = MetricSpec.identity(3).toOption.get
    PreparedPair.from(view(x), view(y), PreprocessSpec.Center, PreprocessSpec.Center, Some(badMetric)) match
      case Left(MultivarError.MetricShapeMismatch(IndexAxis.Row, 4, 3)) => ()
      case other =>
        fail(s"expected metric shape mismatch, got $other")

  test("fromPredictive rejects a response scale that cannot be inverted; from does not"):
    val collapsing = PreprocessSpec.multiplyColumns(Vector(0.0, 2.0)).toOption.get
    PreparedPair.from(view(x), view(y), PreprocessSpec.Center, collapsing) match
      case Right(prepared) =>
        prepared.invertibleResponse match
          case Left(MultivarError.InvalidMap(_)) => ()
          case other => fail(s"ordinary prep must not claim invertibility, got $other")
      case other =>
        fail(s"ordinary from should accept a zero response scale, got $other")

    PreparedPair.fromPredictive(view(x), view(y), PreprocessSpec.Center, collapsing) match
      case Left(MultivarError.NonInvertibleValue("affine inverse scale", 0, 0.0)) => ()
      case other =>
        fail(s"expected zero response scale rejected at predictive fit, got $other")

  test("fromPredictive exposes an invertible response and finite dense working matrices"):
    val prepared = PreparedPair
      .fromPredictive(view(x), view(y), PreprocessSpec.Center, PreprocessSpec.Standardize())
      .toOption
      .get
    assertEquals(prepared.moments.sampleCount, 4)
    assertEquals(prepared.moments.xFeatures, 2)
    assertEquals(prepared.moments.yFeatures, 2)
    assertEqualsDouble(prepared.moments.covarianceScale, 1.0 / 3.0, 0.0)
    assert(prepared.invertibleResponse.isRight)

    val (xd, yd) = prepared.workingDense().toOption.get
    assertEquals(xd.rows, 4)
    assertEquals(yd.cols, 2)
    var col = 0
    while col < xd.cols do
      var sum = 0.0
      var row = 0
      while row < xd.rows do
        sum += xd(row, col)
        row += 1
      assertEqualsDouble(sum, 0.0, 1e-12)
      col += 1

  test("covarianceScale is 1 for a single sample and workingDense rejects non-finite values"):
    val oneX = GaleNumerics.matrixFromRows(Vector(Vector(1.0, 2.0)))
    val oneY = GaleNumerics.matrixFromRows(Vector(Vector(3.0, 4.0)))
    val prepared = PreparedPair.from(view(oneX), view(oneY), PreprocessSpec.Pass, PreprocessSpec.Pass).toOption.get
    assertEqualsDouble(prepared.moments.covarianceScale, 1.0, 0.0)

    val dirty = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, Double.NaN),
        Vector(2.0, 3.0),
        Vector(4.0, 5.0),
        Vector(6.0, 7.0)
      )
    )
    PreparedPair.from(view(dirty), view(y), PreprocessSpec.Pass, PreprocessSpec.Pass).flatMap(_.workingDense()) match
      case Left(MultivarError.NonFiniteValue(_, _, _)) => ()
      case other =>
        fail(s"expected non-finite rejection, got $other")
