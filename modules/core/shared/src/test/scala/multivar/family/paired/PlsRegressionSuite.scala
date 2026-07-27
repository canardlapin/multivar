package multivar
package family.paired

import multivar.core.*
import multivar.workflow.*

import gale.linalg.DMat
import gale.linalg.DVec

/** Laws for PLS regression plus parity against R package `pls`
  * (`plsr(..., method = "simpls")` fixtures).
  */
class PlsRegressionSuite extends munit.FunSuite:

  import PlsRegressionRReferenceFixtures as R

  private def k(value: Int): ComponentCount =
    ComponentCount.unsafe(value)

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

  private def assertVectorClose(actual: DVec, expected: DVec, tol: Double): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), tol)
      i += 1

  private def canonicalizeColumns(weights: DMat, scores: DMat, xLoadings: DMat, yLoadings: DMat)
      : (DMat, DMat, DMat, DMat) =
    val w = weights.copyData
    val s = scores.copyData
    val xl = xLoadings.copyData
    val yl = yLoadings.copyData
    var col = 0
    while col < weights.cols do
      var best = 0
      var row = 1
      while row < weights.rows do
        if Math.abs(w(row * weights.cols + col)) > Math.abs(w(best * weights.cols + col)) then best = row
        row += 1
      if w(best * weights.cols + col) < 0.0 then
        row = 0
        while row < weights.rows do
          w(row * weights.cols + col) = -w(row * weights.cols + col)
          row += 1
        row = 0
        while row < scores.rows do
          s(row * scores.cols + col) = -s(row * scores.cols + col)
          row += 1
        row = 0
        while row < xLoadings.rows do
          xl(row * xLoadings.cols + col) = -xl(row * xLoadings.cols + col)
          row += 1
        row = 0
        while row < yLoadings.rows do
          yl(row * yLoadings.cols + col) = -yl(row * yLoadings.cols + col)
          row += 1
      col += 1
    (
      GaleNumerics.matrixFromRowMajor(weights.rows, weights.cols, w),
      GaleNumerics.matrixFromRowMajor(scores.rows, scores.cols, s),
      GaleNumerics.matrixFromRowMajor(xLoadings.rows, xLoadings.cols, xl),
      GaleNumerics.matrixFromRowMajor(yLoadings.rows, yLoadings.cols, yl)
    )

  private def assertPlsrParity(ref: R.PlsPlsrReference, tol: Double = 1e-9): PlsRegressionFit =
    val fit = PlsRegression.fit(ref.x, ref.y, ref.components).toOption.get
    val (rot, scores, xLoad, yLoad) =
      canonicalizeColumns(fit.xRotations, fit.xScores, fit.xLoadings, fit.yLoadings)

    assertEquals(fit.diagnostics.algorithm, PlsAlgorithm.Simpls)
    assertEquals(fit.effectiveComponents, ref.components)
    assertMatrixClose(fit.coefficients, ref.coefficients, tol)
    assertMatrixClose(fit.workingCoefficients, ref.coefficients, tol)
    assertVectorClose(fit.intercept, ref.intercept, tol)
    assertMatrixClose(rot, ref.xRotations, tol)
    assertMatrixClose(scores, ref.xScores, tol)
    assertMatrixClose(xLoad, ref.xLoadings, tol)
    assertMatrixClose(yLoad, ref.yLoadings, tol)
    assertMatrixClose(fit.predict(ref.x).toOption.get, ref.trainingPredictions, tol)
    assertMatrixClose(fit.predict(ref.newX).toOption.get, ref.newPredictions, tol)
    assertMatrixClose(fit.transform(ref.x).toOption.get, fit.xScores, 1e-12)
    fit

  private def predictByHand(x: DMat, coefficients: DMat, intercept: DVec): DMat =
    val out = Array.ofDim[Double](x.rows * coefficients.cols)
    var row = 0
    while row < x.rows do
      var col = 0
      while col < coefficients.cols do
        var acc = intercept(col)
        var feature = 0
        while feature < x.cols do
          acc += x(row, feature) * coefficients(feature, col)
          feature += 1
        out(row * coefficients.cols + col) = acc
        col += 1
      row += 1
    GaleNumerics.matrixFromRowMajor(x.rows, coefficients.cols, out)

  private def frobenius(left: DMat, right: DMat): Double =
    var acc = 0.0
    var row = 0
    while row < left.rows do
      var col = 0
      while col < left.cols do
        val d = left(row, col) - right(row, col)
        acc += d * d
        col += 1
      row += 1
    Math.sqrt(acc)

  private def columnMeans(matrix: DMat): DVec =
    val out = new Array[Double](matrix.cols)
    var col = 0
    while col < matrix.cols do
      var sum = 0.0
      var row = 0
      while row < matrix.rows do
        sum += matrix(row, col)
        row += 1
      out(col) = sum / matrix.rows.toDouble
      col += 1
    GaleNumerics.vectorFromArray(out)

  private def centerColumns(matrix: DMat): DMat =
    val means = columnMeans(matrix)
    val out = matrix.copyData
    var row = 0
    while row < matrix.rows do
      var col = 0
      while col < matrix.cols do
        out(row * matrix.cols + col) -= means(col)
        col += 1
      row += 1
    GaleNumerics.matrixFromRowMajor(matrix.rows, matrix.cols, out)

  private def solveSymmetric(gram: DMat, rhs: DMat): DMat =
    import gale.linalg.CholeskyOptions
    gram.cholesky(CholeskyOptions()).toOption.get.solve(rhs).toOption.get

  // --- R package pls parity (plsr, method = "simpls") ---

  test("R pls::plsr multivariate SIMPLS matches coefficients, latents, and predictions"):
    assertPlsrParity(R.multivariate)

  test("R pls::plsr univariate multi-component SIMPLS matches"):
    assertPlsrParity(R.univariate)

  test("R pls::plsr wide-response SIMPLS matches (m > p eigen branch)"):
    assertPlsrParity(R.wideResponse)

  test("R pls::plsr yarn NIR subset matches"):
    assertPlsrParity(R.yarnSubset)

  test("R pls training predictions equal x * coefficients + intercept"):
    val ref = R.multivariate
    val fit = assertPlsrParity(ref)
    assertMatrixClose(
      predictByHand(ref.x, fit.coefficients, fit.intercept),
      ref.trainingPredictions,
      1e-9
    )
    assertMatrixClose(
      predictByHand(ref.newX, fit.coefficients, fit.intercept),
      ref.newPredictions,
      1e-9
    )

  // --- Internal laws ---

  test("transform(trainingX) equals stored xScores"):
    val ref = R.multivariate
    val fit = PlsRegression.fit(ref.x, ref.y, ref.components).toOption.get
    assertMatrixClose(fit.transform(ref.x).toOption.get, fit.xScores, 1e-12)

  test("predict(x) equals x * coefficients + intercept in raw coordinates"):
    val ref = R.multivariate
    val fit = PlsRegression.fit(ref.x, ref.y, ref.components).toOption.get
    assertMatrixClose(
      fit.predict(ref.x).toOption.get,
      predictByHand(ref.x, fit.coefficients, fit.intercept),
      1e-10
    )

  test("full-rank SIMPLS predictions converge to OLS"):
    val x = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 0.5),
        Vector(2.0, -1.0),
        Vector(0.0, 1.5),
        Vector(-1.0, 0.0),
        Vector(3.0, 2.0),
        Vector(1.5, -0.5),
        Vector(-0.5, 1.0),
        Vector(2.5, 0.0)
      )
    )
    val y = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.5, 0.0),
        Vector(0.5, 1.0),
        Vector(2.0, -0.5),
        Vector(-1.0, 0.5),
        Vector(4.0, 1.5),
        Vector(0.0, -1.0),
        Vector(1.0, 0.5),
        Vector(2.0, 0.0)
      )
    )
    val fit = PlsRegression.fit(x, y, components = 2).toOption.get
    val xc = centerColumns(x)
    val yc = centerColumns(y)
    val olsWorking = solveSymmetric(GaleNumerics.crossProduct(xc), GaleNumerics.transposeMultiply(xc, yc))
    val mx = columnMeans(x)
    val my = columnMeans(y)
    val intercept = Array.ofDim[Double](y.cols)
    var col = 0
    while col < y.cols do
      var acc = my(col)
      var feature = 0
      while feature < x.cols do
        acc -= mx(feature) * olsWorking(feature, col)
        feature += 1
      intercept(col) = acc
      col += 1
    val olsPred = predictByHand(x, olsWorking, GaleNumerics.vectorFromArray(intercept))
    assert(frobenius(fit.workingCoefficients, olsWorking) < 1e-8)
    assertMatrixClose(fit.predict(x).toOption.get, olsPred, 1e-8)

  test("score columns are orthonormal under SIMPLS"):
    val ref = R.multivariate
    val fit = PlsRegression.fit(ref.x, ref.y, ref.components).toOption.get
    val gram = GaleNumerics.crossProduct(fit.xScores)
    var row = 0
    while row < gram.rows do
      var col = 0
      while col < gram.cols do
        val expected = if row == col then 1.0 else 0.0
        assertEqualsDouble(gram(row, col), expected, 1e-10)
        col += 1
      row += 1

  test("component ordering is deterministic across repeated fits"):
    val ref = R.multivariate
    val first = PlsRegression.fit(ref.x, ref.y, ref.components).toOption.get
    val second = PlsRegression.fit(ref.x, ref.y, ref.components).toOption.get
    assertMatrixClose(first.coefficients, second.coefficients, 0.0)
    assertMatrixClose(first.xRotations, second.xRotations, 0.0)
    assertMatrixClose(first.xScores, second.xScores, 0.0)

  test("unequal rows and over-rank requests are rejected"):
    val x = GaleNumerics.matrixFromRows(Vector(Vector(1.0, 2.0), Vector(3.0, 4.0)))
    val y = GaleNumerics.matrixFromRows(Vector(Vector(1.0), Vector(2.0), Vector(3.0)))
    assert(PlsRegression.fit(x, y, 1).isLeft)
    val okY = GaleNumerics.matrixFromRows(Vector(Vector(1.0), Vector(2.0)))
    // SIMPLS is capped by min(n, p), not by response width.
    assertEquals(
      PlsRegression.fit(x, okY, 3).left.toOption.get,
      MultivarError.InvalidComponentRequest(3, 2)
    )

  test("PairedMultivarEstimator.PlsRegression stays off PairedProgramKind"):
    val estimator = PairedMultivarEstimator.PlsRegression(k(1))
    assertEquals(estimator.label, "pls-simpls")
    assertEquals(estimator.method, PairedLatentMethod.PlsRegression(PlsAlgorithm.Simpls))
    val x = SampleByFeatureInput.of("pls-x", 12, 4).toOption.get
    val y = SampleByFeatureInput.of("pls-y", 12, 2).toOption.get
    val paired = PairedSampleByFeatureInput.of(x, y).toOption.get
    val plan = PairedMultivarPlan.of("pls-pair", paired, estimator).toOption.get
    assertEquals(plan.estimator.label, "pls-simpls")

  test("PairedMultivarPlan allows univariate multi-component PLS"):
    val x = SampleByFeatureInput.of("pls-uni-x", 12, 4).toOption.get
    val y = SampleByFeatureInput.of("pls-uni-y", 12, 1).toOption.get
    val paired = PairedSampleByFeatureInput.of(x, y).toOption.get
    val plan = PairedMultivarPlan.of("pls-uni", paired, PairedMultivarEstimator.PlsRegression(k(2)))
    assert(plan.isRight, s"expected univariate multi-component PLS plan, got $plan")

  test("raw coefficients satisfy the original-coordinate prediction law under invertible affine preprocessing"):
    val schemes = Vector(
      (PreprocessSpec.Center, PreprocessSpec.Center),
      (PreprocessSpec.Standardize(), PreprocessSpec.Standardize()),
      (
        PreprocessSpec.multiplyColumns(Vector(2.0, 0.5, 1.5, 1.0)).toOption.get,
        PreprocessSpec.multiplyColumns(Vector(3.0, 0.25)).toOption.get
      )
    )
    val ref = R.multivariate
    schemes.foreach { case (xPre, yPre) =>
      val fit = PlsRegression
        .fit(ref.x, ref.y, ref.components, xPre, yPre)
        .toOption
        .get
      assertMatrixClose(
        fit.predict(ref.x).toOption.get,
        predictByHand(ref.x, fit.coefficients, fit.intercept),
        1e-9
      )
      assertMatrixClose(
        fit.predict(ref.newX).toOption.get,
        predictByHand(ref.newX, fit.coefficients, fit.intercept),
        1e-9
      )
    }

  test("predictWorking and transform agree with the working-space maps under Standardize"):
    val ref = R.multivariate
    val xPre = PreprocessSpec.Standardize()
    val yPre = PreprocessSpec.Standardize()
    val fit = PlsRegression.fit(ref.x, ref.y, ref.components, xPre, yPre).toOption.get
    val fittedX = xPre.fit(MatrixView.dense(ref.x)).toOption.get
    val xWorking = fittedX.transform(MatrixView.dense(ref.x)).toOption.get.toDense().toOption.get
    val scoresByHand = GaleNumerics.multiply(xWorking, fit.xRotations)
    assertMatrixClose(fit.transform(ref.x).toOption.get, scoresByHand, 1e-9)
    val workingPred = GaleNumerics.multiply(xWorking, fit.workingCoefficients)
    assertMatrixClose(fit.predictWorking(ref.x).toOption.get, workingPred, 1e-9)
    assertMatrixClose(
      fit.predict(ref.x).toOption.get,
      predictByHand(ref.x, fit.coefficients, fit.intercept),
      1e-9
    )

  test("PLS rejects a response scale it could not undo when predicting"):
    val ref = R.multivariate
    val collapsing = PreprocessSpec.multiplyColumns(Vector.fill(ref.y.cols)(0.0)).toOption.get
    PlsRegression.fit(ref.x, ref.y, ref.components, PreprocessSpec.Center, collapsing) match
      case Left(MultivarError.NonInvertibleValue("affine inverse scale", 0, 0.0)) => ()
      case other => fail(s"expected zero response scale rejected at fit time, got $other")
