package multivar
package family.paired

import multivar.core.*
import multivar.family.paired.*
import multivar.workflow.*
import multivar.advanced.*

import gale.linalg.DMat
import gale.linalg.DVec

class PairedLatentSuite extends munit.FunSuite:

  import PairedLatentRReferenceFixtures as R

  private def plscOp(fit: PlscFit) = PlscFit.operatorOf(fit)
  private def ccaOp(fit: CcaFit) = CcaFit.operatorOf(fit)
  private def rrrOp(fit: ReducedRankRegressionFit) = ReducedRankRegressionFit.operatorOf(fit)

  private final case class CanonicalPaired(
      xWeights: DMat,
      yWeights: DMat,
      xScores: DMat,
      yScores: DMat
  )

  private def k(value: Int): ComponentCount =
    ComponentCount.unsafe(value)

  private def assertMatrixClose(actual: DMat, expected: Vector[Vector[Double]], tol: Double): Unit =
    assertEquals(actual.rows, expected.length)
    assertEquals(actual.cols, expected.headOption.map(_.length).getOrElse(0))
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assertEqualsDouble(actual(row, col), expected(row)(col), tol)
        col += 1
      row += 1

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

  private def canonicalize(
      xWeights: DMat,
      yWeights: DMat,
      xScores: DMat,
      yScores: DMat
  ): CanonicalPaired =
    val xWeightData = xWeights.copyData
    val yWeightData = yWeights.copyData
    val xScoreData = xScores.copyData
    val yScoreData = yScores.copyData
    var col = 0
    while col < xWeights.cols do
      var best = 0
      var row = 1
      while row < xWeights.rows do
        if Math.abs(xWeightData(row * xWeights.cols + col)) > Math.abs(xWeightData(best * xWeights.cols + col)) then best = row
        row += 1
      if xWeightData(best * xWeights.cols + col) < 0.0 then
        row = 0
        while row < xWeights.rows do
          xWeightData(row * xWeights.cols + col) = -xWeightData(row * xWeights.cols + col)
          row += 1
        row = 0
        while row < yWeights.rows do
          yWeightData(row * yWeights.cols + col) = -yWeightData(row * yWeights.cols + col)
          row += 1
        row = 0
        while row < xScores.rows do
          xScoreData(row * xScores.cols + col) = -xScoreData(row * xScores.cols + col)
          row += 1
        row = 0
        while row < yScores.rows do
          yScoreData(row * yScores.cols + col) = -yScoreData(row * yScores.cols + col)
          row += 1
      col += 1
    CanonicalPaired(
      GaleNumerics.matrixFromRowMajor(xWeights.rows, xWeights.cols, xWeightData),
      GaleNumerics.matrixFromRowMajor(yWeights.rows, yWeights.cols, yWeightData),
      GaleNumerics.matrixFromRowMajor(xScores.rows, xScores.cols, xScoreData),
      GaleNumerics.matrixFromRowMajor(yScores.rows, yScores.cols, yScoreData)
    )

  private def assertPairedReference(actual: CanonicalPaired, expected: R.PairedReference, tol: Double): Unit =
    assertMatrixClose(actual.xWeights, expected.xWeights, tol)
    assertMatrixClose(actual.yWeights, expected.yWeights, tol)
    assertMatrixClose(actual.xScores, expected.xScores, tol)
    assertMatrixClose(actual.yScores, expected.yScores, tol)

  test("every representable paired estimator is executable") {
    val x = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(1.0, 0.0),
          Vector(0.0, 2.0),
          Vector(3.0, 1.0),
          Vector(4.0, -1.0),
          Vector(-2.0, 1.5)
        )
      )
    )
    val y = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(2.0, 1.0),
          Vector(1.0, -1.0),
          Vector(0.0, 3.0),
          Vector(5.0, 2.0),
          Vector(-1.0, 0.5)
        )
      )
    )
    val estimators = Vector[PairedMultivarEstimator](
      PairedMultivarEstimator.Plsc(k(1)),
      PairedMultivarEstimator.Cca(k(1)),
      PairedMultivarEstimator.ReducedRankRegression(k(1)),
      PairedMultivarEstimator.PlsRegression(k(1))
    )

    def fitKind(estimator: PairedMultivarEstimator): Either[MultivarError, String] =
      estimator match
        case PairedMultivarEstimator.Plsc(components, xPreprocessing, yPreprocessing) =>
          Plsc.fit(x, y, components, xPreprocessing, yPreprocessing).map(fit => plscOp(fit).diagnostics.kind.label)
        case PairedMultivarEstimator.Cca(components, regularization, xPreprocessing, yPreprocessing) =>
          Cca.fitRegularized(x, y, components, regularization, xPreprocessing, yPreprocessing)
            .map(fit => ccaOp(fit).diagnostics.kind.label)
        case PairedMultivarEstimator.ReducedRankRegression(
              components,
              regularization,
              direction,
              xPreprocessing,
              yPreprocessing
            ) =>
          ReducedRankRegression
            .fit(x, y, components, regularization, direction, xPreprocessing, yPreprocessing)
            .map(fit => rrrOp(fit).diagnostics.kind.label)
        case PairedMultivarEstimator.PlsRegression(components, algorithm, xPreprocessing, yPreprocessing) =>
          PlsRegression
            .fit(x, y, components, xPreprocessing, yPreprocessing, algorithm)
            .map(fit => fit.diagnostics.algorithm.label)
            .map(label => s"pls-$label")

    estimators.foreach { estimator =>
      val actual = fitKind(estimator)
      assert(actual.isRight, s"${estimator.label} should be executable, got $actual")
      assertEquals(actual.toOption.get, estimator.label)
    }
  }

  test("PLSC matches the R cross-covariance SVD fixture") {
    val ref = R.plsc
    val fit = Plsc.fit(MatrixView.dense(ref.x), MatrixView.dense(ref.y), k(ref.components)).toOption.get
    val actual = canonicalize(
      plscOp(fit).sourceWeights.toOption.get,
      plscOp(fit).targetWeights.toOption.get,
      fit.xScores,
      fit.yScores
    )

    assertVectorClose(fit.svdResult.singularValues, ref.singularValues, 1e-9)
    assertEquals(plscOp(fit).diagnostics.kind, PairedProgramKind.Plsc)
    assertPairedReference(actual, ref, 1e-9)
  }

  test("regularized CCA matches the R ridge whitening fixture and base cancor correlations") {
    val ref = R.cca
    val regularization = CcaRegularization.asymmetric(ref.xRidge, ref.yRidge).toOption.get
    val fit = Cca.fitRegularized(MatrixView.dense(ref.x), MatrixView.dense(ref.y), k(ref.components), regularization).toOption.get
    val actual = canonicalize(
      ccaOp(fit).sourceWeights.toOption.get,
      ccaOp(fit).targetWeights.toOption.get,
      fit.xScores,
      fit.yScores
    )

    assertVectorClose(fit.svdResult.singularValues, ref.singularValues, 1e-9)
    ccaOp(fit).diagnostics.kind match
      case PairedProgramKind.Cca(value) => assertEquals(value, regularization)
      case other                        => fail(s"expected CCA kind, got $other")
    assertPairedReference(actual, ref, 1e-9)

    val zeroRidge = CcaRegularization.asymmetric(0.0, 0.0).toOption.get
    val baseFit = Cca.fitRegularized(MatrixView.dense(ref.x), MatrixView.dense(ref.y), k(ref.components), zeroRidge).toOption.get
    assertVectorClose(baseFit.svdResult.singularValues, ref.baseCancorCorrelations.get, 1e-9)
  }

  test("PLSC with Standardize matches Pass on hand-standardized inputs up to sign"):
    val x = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 2.0, 0.5),
        Vector(2.0, 1.0, -0.5),
        Vector(0.0, 3.0, 1.0),
        Vector(-1.0, 0.5, 0.0),
        Vector(3.0, -1.0, 2.0),
        Vector(1.5, 0.0, -1.0)
      )
    )
    val y = GaleNumerics.matrixFromRows(
      Vector(
        Vector(2.0, 1.0),
        Vector(1.0, -1.0),
        Vector(0.0, 3.0),
        Vector(5.0, 2.0),
        Vector(-1.0, 0.5),
        Vector(0.5, 1.5)
      )
    )
    val direct = Plsc.fit(x, y, 1, PreprocessSpec.Standardize(), PreprocessSpec.Standardize()).toOption.get
    val xPre = PreprocessSpec.Standardize().fit(MatrixView.dense(x)).toOption.get
    val yPre = PreprocessSpec.Standardize().fit(MatrixView.dense(y)).toOption.get
    val xHat = xPre.transform(MatrixView.dense(x)).toOption.get.toDense().toOption.get
    val yHat = yPre.transform(MatrixView.dense(y)).toOption.get.toDense().toOption.get
    val viaPass = Plsc.fit(xHat, yHat, 1, PreprocessSpec.Pass, PreprocessSpec.Pass).toOption.get
    val left = canonicalize(direct.xWeights, direct.yWeights, direct.xScores, direct.yScores)
    val right = canonicalize(viaPass.xWeights, viaPass.yWeights, viaPass.xScores, viaPass.yScores)
    assertMatrixClose(left.xWeights, right.xWeights, 1e-9)
    assertMatrixClose(left.yWeights, right.yWeights, 1e-9)
    assertMatrixClose(left.xScores, right.xScores, 1e-9)
    assertMatrixClose(left.yScores, right.yScores, 1e-9)

  test("PLSC uses the shared row metric in the paired cross operator") {
    val x = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 0.0),
        Vector(0.0, 2.0),
        Vector(3.0, 1.0),
        Vector(4.0, -1.0),
        Vector(-2.0, 1.5)
      )
    )
    val y = GaleNumerics.matrixFromRows(
      Vector(
        Vector(2.0, 1.0),
        Vector(1.0, -1.0),
        Vector(0.0, 3.0),
        Vector(5.0, 2.0),
        Vector(-1.0, 0.5)
      )
    )
    val weights = Vector(1.0, 2.0, 0.5, 1.5, 0.75)
    val rowMetric = MetricSpec.diagonal(DVec.fromSeq(weights)).toOption.get
    val xView = MatrixView.dense(x)
    val yView = MatrixView.dense(y)
    val fit = Plsc.fit(xView, yView, k(2), rowMetric = Some(rowMetric)).toOption.get
    val unweighted = Plsc.fit(xView, yView, k(2)).toOption.get
    val xp = PreprocessSpec.Center.fit(xView).toOption.get.transform(xView).toOption.get
    val yp = PreprocessSpec.Center.fit(yView).toOption.get.transform(yView).toOption.get
    val cross = MatrixOps.scale(
      weightedCrossProduct(xp.toDense().toOption.get, yp.toDense().toOption.get, weights),
      1.0 / (x.rows - 1)
    )
    val expected = DenseSolvers.svd.decompose(MatrixView.dense(cross), k(2)).toOption.get

    assertVectorClose(fit.svdResult.singularValues, expected.singularValues, 1e-9)
    assert(Math.abs(fit.svdResult.singularValues(0) - unweighted.svdResult.singularValues(0)) > 1e-3)
  }

  test("CCA uses the shared row metric in covariance and cross-covariance operators") {
    val x = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 0.0),
        Vector(0.0, 2.0),
        Vector(3.0, 1.0),
        Vector(4.0, -1.0),
        Vector(-2.0, 1.5)
      )
    )
    val y = GaleNumerics.matrixFromRows(
      Vector(
        Vector(2.0, 1.0),
        Vector(1.0, -1.0),
        Vector(0.0, 3.0),
        Vector(5.0, 2.0),
        Vector(-1.0, 0.5)
      )
    )
    val ridge = 0.2
    val weights = Vector(1.0, 2.0, 0.5, 1.5, 0.75)
    val rowMetric = MetricSpec.diagonal(DVec.fromSeq(weights)).toOption.get
    val xView = MatrixView.dense(x)
    val yView = MatrixView.dense(y)
    val fit = Cca.fit(xView, yView, k(2), ridge = ridge, rowMetric = Some(rowMetric)).toOption.get
    val unweighted = Cca.fit(xView, yView, k(2), ridge = ridge).toOption.get
    val xp = PreprocessSpec.Center.fit(xView).toOption.get.transform(xView).toOption.get
    val yp = PreprocessSpec.Center.fit(yView).toOption.get.transform(yView).toOption.get
    val scale = 1.0 / (x.rows - 1)
    val xpDense = xp.toDense().toOption.get
    val ypDense = yp.toDense().toOption.get
    val cxx = MatrixOps.addRidge(MatrixOps.scale(weightedCrossProduct(xpDense, xpDense, weights), scale), ridge)
    val cyy = MatrixOps.addRidge(MatrixOps.scale(weightedCrossProduct(ypDense, ypDense, weights), scale), ridge)
    val cxy = MatrixOps.scale(weightedCrossProduct(xpDense, ypDense, weights), scale)
    val wx = MatrixOps.inverseSquareRoot(cxx, DenseSolvers.symmetricEigen, 1e-12).toOption.get
    val wy = MatrixOps.inverseSquareRoot(cyy, DenseSolvers.symmetricEigen, 1e-12).toOption.get
    val expectedOperator = GaleNumerics.multiply(GaleNumerics.multiply(wx, cxy), wy)
    val expected = DenseSolvers.svd.decompose(MatrixView.dense(expectedOperator), k(2)).toOption.get

    assertVectorClose(fit.svdResult.singularValues, expected.singularValues, 1e-9)
    assert(Math.abs(fit.svdResult.singularValues(1) - unweighted.svdResult.singularValues(1)) > 1e-3)
  }

  test("weighted PLSC matches an externally computed weighted cross-product SVD") {
    val x = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 0.0),
        Vector(0.0, 2.0),
        Vector(3.0, 1.0),
        Vector(4.0, -1.0),
        Vector(-2.0, 1.5)
      )
    )
    val y = GaleNumerics.matrixFromRows(
      Vector(
        Vector(2.0, 1.0),
        Vector(1.0, -1.0),
        Vector(0.0, 3.0),
        Vector(5.0, 2.0),
        Vector(-1.0, 0.5)
      )
    )
    val weights = Vector(1.0, 2.0, 0.5, 1.5, 0.75)
    val rowMetric = MetricSpec.diagonal(DVec.fromSeq(weights)).toOption.get
    val fit = Plsc.fit(MatrixView.dense(x), MatrixView.dense(y), k(2), rowMetric = Some(rowMetric)).toOption.get

    // External reference: column-center each table, then build X' D Y / (n - 1)
    // entrywise with plain loops — no DualityKernels, no metric machinery.
    val n = x.rows
    def centered(data: DMat): DMat =
      val out = data.copyData
      var col = 0
      while col < data.cols do
        var mean = 0.0
        var row = 0
        while row < n do
          mean += data(row, col)
          row += 1
        mean /= n
        row = 0
        while row < n do
          out(row * data.cols + col) -= mean
          row += 1
        col += 1
      GaleNumerics.matrixFromRowMajor(data.rows, data.cols, out)
    val xc = centered(x)
    val yc = centered(y)
    val crossData = new Array[Double](x.cols * y.cols)
    var p = 0
    while p < x.cols do
      var q = 0
      while q < y.cols do
        var acc = 0.0
        var row = 0
        while row < n do
          acc += weights(row) * xc(row, p) * yc(row, q)
          row += 1
        crossData(p * y.cols + q) = acc / (n - 1)
        q += 1
      p += 1
    val cross = GaleNumerics.matrixFromRowMajor(x.cols, y.cols, crossData)
    val expected = DenseSolvers.svd.decompose(MatrixView.dense(cross), k(2)).toOption.get

    assertVectorClose(fit.svdResult.singularValues, expected.singularValues, 1e-9)
    val actual = canonicalize(
      plscOp(fit).sourceWeights.toOption.get,
      plscOp(fit).targetWeights.toOption.get,
      fit.xScores,
      fit.yScores
    )
    val expectedScores = canonicalize(
      expected.u,
      expected.v,
      GaleNumerics.multiply(xc, expected.u),
      GaleNumerics.multiply(yc, expected.v)
    )
    assertMatrixClose(actual.xWeights, expectedScores.xWeights, 1e-9)
    assertMatrixClose(actual.yWeights, expectedScores.yWeights, 1e-9)
    assertMatrixClose(actual.xScores, expectedScores.xScores, 1e-9)
    assertMatrixClose(actual.yScores, expectedScores.yScores, 1e-9)
  }

  test("RRR x-to-y matches the R reduced-rank regression fixture") {
    val ref = R.rrr
    val fit = ReducedRankRegression.fit(MatrixView.dense(ref.x), MatrixView.dense(ref.y), k(ref.components)).toOption.get
    val source = ReducedRankRegressionFit.sourceOf(fit)
    val target = ReducedRankRegressionFit.targetOf(fit)
    val actual = canonicalize(
      source.frame.weights.toDense.toOption.get,
      target.frame.weights.toDense.toOption.get,
      source.trainingValues,
      target.trainingValues
    )

    assertVectorClose(ReducedRankRegressionFit.operatorOf(fit).result.singularValues, ref.singularValues, 1e-9)
    assertPairedReference(actual, ref, 1e-9)
    assertMatrixClose(
      ReducedRankRegressionFit.unconstrainedWorkingCoefficients(fit),
      ref.fullCoefficient.get,
      1e-9
    )
    assertMatrixClose(fit.workingCoefficients, ref.workingCoefficient.get, 1e-9)
    assertMatrixClose(fit.predictWorking(MatrixView.dense(ref.x)).toOption.get, ref.predictedWorking.get, 1e-9)
    assertMatrixClose(fit.predict(MatrixView.dense(ref.x)).toOption.get, ref.predicted.get, 1e-9)
  }

  test("RRR raw coefficients satisfy the original-coordinate prediction law under invertible affine preprocessing") {
    val x = GaleNumerics.matrixFromRows(
      Vector(
        Vector(2.5, 2.4, 0.5),
        Vector(0.5, 0.7, -0.1),
        Vector(2.2, 2.9, 0.8),
        Vector(1.9, 2.2, 0.3),
        Vector(3.1, 3.0, 1.0),
        Vector(1.2, 1.1, 0.0)
      )
    )
    val y = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 0.5),
        Vector(0.2, -0.3),
        Vector(1.1, 0.7),
        Vector(0.8, 0.4),
        Vector(1.4, 0.9),
        Vector(0.3, -0.1)
      )
    )
    val schemes = Vector(
      (PreprocessSpec.Center, PreprocessSpec.Center),
      (PreprocessSpec.Standardize(), PreprocessSpec.Standardize()),
      (
        PreprocessSpec.multiplyColumns(Vector(2.0, 0.5, 1.5)).toOption.get,
        PreprocessSpec.multiplyColumns(Vector(3.0, 0.25)).toOption.get
      )
    )

    schemes.foreach { case (xPre, yPre) =>
      val fit = ReducedRankRegression
        .fit(x, y, components = 2, RegressionRegularization.Ols, xPre, yPre)
        .toOption
        .get
      val predicted = fit.predict(x).toOption.get
      val fromRaw = applyRaw(x, fit.coefficients, fit.intercept)
      assertMatrixClose(predicted, fromRaw, 1e-9)
    }
  }

  private def applyRaw(x: DMat, coefficients: DMat, intercept: DVec): Vector[Vector[Double]] =
    (0 until x.rows).map { row =>
      (0 until coefficients.cols).map { col =>
        var value = intercept(col)
        var feature = 0
        while feature < coefficients.rows do
          value += x(row, feature) * coefficients(feature, col)
          feature += 1
        value
      }.toVector
    }.toVector

  test("RRR rejects a response scale it could not undo when predicting") {
    val ref = R.rrr
    val collapsing =
      PreprocessSpec.multiplyColumns(Vector.fill(ref.y.cols)(0.0)).toOption.get
    val usable =
      PreprocessSpec.multiplyColumns(Vector.fill(ref.y.cols)(2.0)).toOption.get

    val rejected = ReducedRankRegression.fit(
      ref.x,
      ref.y,
      ref.components,
      RegressionRegularization.Ols,
      PreprocessSpec.Center,
      collapsing
    )
    rejected match
      case Left(MultivarError.NonInvertibleValue("affine inverse scale", 0, 0.0)) => ()
      case other => fail(s"expected the zero response scale to be rejected at fit time, got $other")

    val accepted = ReducedRankRegression.fit(
      ref.x,
      ref.y,
      ref.components,
      RegressionRegularization.Ols,
      PreprocessSpec.Center,
      usable
    )
    assert(accepted.flatMap(_.predict(MatrixView.dense(ref.x))).isRight)
  }

  test("CCA regularization normalizes raw doubles through linalg ridge") {
    val symmetric = CcaRegularization.symmetric(0.25).toOption.get
    val asymmetric = CcaRegularization.asymmetric(0.1, 0.2).toOption.get

    assertEqualsDouble(symmetric.x.value, 0.25, 0.0)
    assertEqualsDouble(symmetric.y.value, 0.25, 0.0)
    assertEqualsDouble(asymmetric.x.value, 0.1, 0.0)
    assertEqualsDouble(asymmetric.y.value, 0.2, 0.0)
    assert(CcaRegularization.symmetric(Double.NaN).isLeft)
    assert(CcaRegularization.asymmetric(0.0, -1.0).isLeft)
  }

  test("regression method records direction and regularization explicitly") {
    val ridge = RegressionRegularization.ridge(0.5).toOption.get
    val method = PairedLatentMethod.ReducedRankRegression(RegressionDirection.XToY, ridge)

    assertEquals(RegressionDirection.XToY.label, "x-to-y")
    assertEquals(ridge.label, "ridge")
    assertEquals(method.label, "rrr")
    assert(RegressionRegularization.ridge(Double.PositiveInfinity).isLeft)
  }

  private def weightedCrossProduct(left: DMat, right: DMat, weights: Vector[Double]): DMat =
    require(left.rows == right.rows && left.rows == weights.length)
    val out = new Array[Double](left.cols * right.cols)
    var source = 0
    while source < left.cols do
      var target = 0
      while target < right.cols do
        var value = 0.0
        var row = 0
        while row < left.rows do
          value += left(row, source) * weights(row) * right(row, target)
          row += 1
        out(source * right.cols + target) = value
        target += 1
      source += 1
    GaleNumerics.matrixFromRowMajor(left.cols, right.cols, out)
