package multivar
package family.paired

import multivar.analysis.*
import multivar.core.*

import gale.linalg.Matrix

class TaskComponentsSuite extends munit.FunSuite:

  private val tol = 1e-8

  private val x = Matrix(6, 2)(
    1.0, 0.0,
    2.0, 1.0,
    3.0, 1.0,
    4.0, 2.0,
    5.0, 3.0,
    6.0, 4.0
  )

  private val y = Matrix(6, 2)(
    2.0, 1.0,
    4.0, 1.5,
    6.0, 2.0,
    8.0, 3.0,
   10.0, 4.0,
   12.0, 5.0
  )

  test("Whole Focus agrees with ReducedRankRegression on predictions and coefficients"):
    val task = TaskComponents.fit(x, y, components = 1).toOption.get
    val rrr = ReducedRankRegression.fit(x, y, components = 1).toOption.get
    assertMatrixClose(task.coefficients, rrr.coefficients, tol)
    assertVectorClose(task.intercept, rrr.intercept, tol)
    assertMatrixClose(task.predict(x).toOption.get, rrr.predict(x).toOption.get, tol)
    assertEquals(task.remainder.rank, 0)
    assertEqualsDouble(task.remainder.value.availableGain, 0.0, tol)

  test("available gain splits across focus and remainder"):
    val groups = Matrix(2, 1)(1.0, 0.0)
    val fit =
      TaskComponents
        .fit(
          x,
          y,
          components = 1,
          resolution = TaskResolution.Target.weights(groups)
        )
        .toOption
        .get
    assertEqualsDouble(
      fit.value.availableGain,
      fit.focus.value.availableGain + fit.remainder.value.availableGain,
      tol
    )
    assertEqualsDouble(
      fit.value.fittedLoss + fit.value.retainedGain,
      fit.value.baselineLoss,
      tol
    )

  test("predict equals baseline plus focus and remainder contributions"):
    val groups = Matrix(2, 2)(
      1.0, 0.0,
      0.0, 1.0
    )
    val fit =
      TaskComponents
        .fit(
          predictors = x,
          responses = y,
          rank = RankBudget.Split(focus = 1, remainder = 1),
          resolution = TaskResolution.Target.structure(groups)
        )
        .toOption
        .get
    val predicted = fit.predict(x).toOption.get
    val base = fit.baseline(x.rows).toOption.get
    val focus = fit.focusContribution(x).toOption.get
    val rem = fit.remainderContribution(x).toOption.get
    var row = 0
    while row < predicted.rows do
      var col = 0
      while col < predicted.cols do
        assertEqualsDouble(
          predicted(row, col),
          base(row, col) + focus(row, col) + rem(row, col),
          tol
        )
        col += 1
      row += 1

  test("Target.weights design equivariance under nonsingular reparameterization"):
    val h = Matrix(2, 1)(1.0, 2.0)
    val ha = Matrix(2, 1)(2.0, 4.0)
    val left =
      TaskComponents.fit(x, y, 1, TaskResolution.Target.weights(h)).toOption.get
    val right =
      TaskComponents.fit(x, y, 1, TaskResolution.Target.weights(ha)).toOption.get
    assertEqualsDouble(left.focus.value.availableGain, right.focus.value.availableGain, tol)
    assertEqualsDouble(left.remainder.value.availableGain, right.remainder.value.availableGain, tol)

  test("Phase 1 rejects Frobenius ridge with nonidentity target metric"):
    val metric = MetricSpec.diagonal(gale.linalg.DVec.fromSeq(Seq(1.0, 2.0))).toOption.get
    val ridge = RegressionRegularization.ridge(1e-3).toOption.get
    val result =
      TaskComponents.fit(
        x,
        y,
        RankBudget.Focus(1),
        TaskResolution.Whole,
        TaskComponents.Options(
          regularization = ridge,
          targetMetric = Some(metric)
        )
      )
    assert(result.isLeft)
    result.left.foreach:
      case MultivarError.InvalidRegularization(_, _, requirement) =>
        assert(requirement.contains("Sylvester"))
      case other =>
        fail(s"expected InvalidRegularization, got $other")

  test("DualSpectrum.Ridge is distinct in provenance from MoorePenrose"):
    val h = Matrix(2, 1)(1.0, 0.0)
    val mp =
      TaskComponents
        .fit(x, y, 1, TaskResolution.Target.structure(h, DualSpectrum.moorePenrose))
        .toOption
        .get
    val ridge =
      TaskComponents
        .fit(x, y, 1, TaskResolution.Target.structure(h, DualSpectrum.Ridge(1e-2)))
        .toOption
        .get
    assertEquals(mp.dualSpectrum.map(_.label), Some("moore-penrose"))
    assertEquals(ridge.dualSpectrum.map(_.label), Some("ridge"))

  test("DualSpectrum.Truncated rejects a rank beyond the task spectrum"):
    val h = Matrix(2, 1)(1.0, 0.0)
    val result =
      TaskComponents.fit(
        x,
        y,
        1,
        TaskResolution.Target.structure(h, DualSpectrum.truncated(components = 3))
      )
    result match
      case Left(MultivarError.InvalidComponentRequest(requested, available)) =>
        assertEquals(requested, 3)
        assertEquals(available, 2)
      case other =>
        fail(s"expected InvalidComponentRequest, got $other")

  test("Source.covariance rejects an observation-axis design with an axis-aware error"):
    val observationDesign = Matrix(6, 1)(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    val result = TaskComponents.fit(x, y, 1, TaskResolution.Source.covariance(observationDesign))
    assert(result.isLeft)
    result.left.foreach:
      case MultivarError.MatrixShapeMismatch(detail) =>
        assert(detail.contains("Source.covariance"))
        assert(detail.contains("2 predictor features"))
        assert(detail.contains("6 rows"))
      case other =>
        fail(s"expected MatrixShapeMismatch, got $other")

  test("Source.covariance and Source.regression are distinct resolvents"):
    val g = Matrix(2, 1)(1.0, 0.0)
    val cov =
      TaskComponents.fit(x, y, 1, TaskResolution.Source.covariance(g)).toOption.get
    val reg =
      TaskComponents.fit(x, y, 1, TaskResolution.Source.regression(g)).toOption.get
    assertEquals(cov.resolution.label, "source covariance")
    assertEquals(reg.resolution.label, "source regression")
    assert(cov.taskSupportRank > 0)
    assertEqualsDouble(
      cov.value.availableGain,
      cov.focus.value.availableGain + cov.remainder.value.availableGain,
      tol
    )
    assertEqualsDouble(
      reg.value.availableGain,
      reg.focus.value.availableGain + reg.remainder.value.availableGain,
      tol
    )
    // Regression retains only the task-supported part of G*, so its focus gain
    // cannot exceed the covariance focus gain for the same generators.
    assert(reg.focus.value.availableGain <= cov.focus.value.availableGain + tol)

  test("Source.regression design equivariance under nonsingular reparameterization"):
    val g = Matrix(2, 1)(1.0, 2.0)
    val ga = Matrix(2, 1)(2.0, 4.0)
    val left =
      TaskComponents.fit(x, y, 1, TaskResolution.Source.regression(g)).toOption.get
    val right =
      TaskComponents.fit(x, y, 1, TaskResolution.Source.regression(ga)).toOption.get
    assertEqualsDouble(left.focus.value.availableGain, right.focus.value.availableGain, tol)
    assertEqualsDouble(left.remainder.value.availableGain, right.remainder.value.availableGain, tol)

  private def assertMatrixClose(left: gale.linalg.DMat, right: gale.linalg.DMat, tolerance: Double): Unit =
    assertEquals(left.rows, right.rows)
    assertEquals(left.cols, right.cols)
    var row = 0
    while row < left.rows do
      var col = 0
      while col < left.cols do
        assertEqualsDouble(left(row, col), right(row, col), tolerance)
        col += 1
      row += 1

  private def assertVectorClose(left: gale.linalg.DVec, right: gale.linalg.DVec, tolerance: Double): Unit =
    assertEquals(left.length, right.length)
    var i = 0
    while i < left.length do
      assertEqualsDouble(left(i), right(i), tolerance)
      i += 1
