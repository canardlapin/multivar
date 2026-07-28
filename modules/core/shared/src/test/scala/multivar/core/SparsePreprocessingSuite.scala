package multivar
package core

import multivar.core.*

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.MutableDVec

class SparsePreprocessingSuite extends munit.FunSuite:

  private val denseRows: Vector[Vector[Double]] =
    Vector(
      Vector(1.0, 0.0, 2.0),
      Vector(0.0, 3.0, 0.0),
      Vector(4.0, 0.0, 5.0),
      Vector(0.0, 6.0, 0.0)
    )

  private def denseMatrix: DMat =
    GaleNumerics.matrixFromRows(denseRows)

  private def sparseView: MatrixView =
    SparseMatrixView.fromRows(denseRows).toOption.get

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

  private def assertVectorClose(actual: DVec, expected: DVec, tol: Double): Unit =
    assertEquals(actual.length, expected.length)
    var index = 0
    while index < actual.length do
      assertEqualsDouble(actual(index), expected(index), tol)
      index += 1

  private def summaryOf(fitted: FittedInvertiblePreprocessor): ColumnAffineSummary =
    fitted match
      case affine: InvertibleColumnAffine => affine.summary
      case other                          => fail(s"expected InvertibleColumnAffine, got $other")

  private def denseTransform(scale: Vector[Double], shift: Vector[Double]): Vector[Vector[Double]] =
    denseRows.map { row =>
      row.zipWithIndex.map { case (value, col) => value * scale(col) + shift(col) }
    }

  test("pass and column scaling preserve sparse storage") {
    val pass = PreprocessSpec.Pass.fit(sparseView).toOption.get
    val passOut = pass.transform(sparseView).toOption.get
    assertEquals(passOut.storage, StorageKind.Sparse)

    val scaleSpec = PreprocessSpec.multiplyColumns(Vector(2.0, 0.5, -1.0)).toOption.get
    val scale = scaleSpec.fit(sparseView).toOption.get
    val scaled = scale.transform(sparseView).toOption.get

    assertEquals(scaled.storage, StorageKind.Sparse)
    assertMatrixClose(
      scaled.toDense(StoragePolicy.AllowDense).toOption.get,
      denseTransform(Vector(2.0, 0.5, -1.0), Vector(0.0, 0.0, 0.0)),
      1e-12
    )
  }

  test("sparse centering is lazy by default and explicit about densification policy") {
    val center = PreprocessSpec.Center.fit(sparseView).toOption.get
    val centered = center.transform(sparseView).toOption.get

    assertEquals(centered.storage, StorageKind.LazyAffine)
    assert(center.transform(sparseView, policy = StoragePolicy.PreserveSparse).isLeft)

    val means = Vector(1.25, 2.25, 1.75)
    assertMatrixClose(
      centered.toDense(StoragePolicy.AllowDense).toOption.get,
      denseTransform(Vector(1.0, 1.0, 1.0), means.map(-_)),
      1e-12
    )

    val denseCentered = center.transform(sparseView, policy = StoragePolicy.AllowDense).toOption.get
    assertEquals(denseCentered.storage, StorageKind.Dense)
  }

  test("standardizing sparse matrices stays lazy and matches dense arithmetic") {
    val fitted = PreprocessSpec.Standardize().fit(sparseView).toOption.get
    val standardized = fitted.transform(sparseView).toOption.get
    val means = Vector(1.25, 2.25, 1.75)
    val sds = Vector(1.8929694486000912, 2.8722813232690143, 2.362907813126304)
    val expected = denseRows.map { row =>
      row.zipWithIndex.map { case (value, col) => (value - means(col)) / sds(col) }
    }

    assertEquals(standardized.storage, StorageKind.LazyAffine)
    assertMatrixClose(standardized.toDense(StoragePolicy.AllowDense).toOption.get, expected, 1e-12)
  }

  test("affine MatrixView algebra matches dense arithmetic") {
    val fitted = PreprocessSpec.Center.fit(sparseView).toOption.get
    val centered = fitted.transform(sparseView).toOption.get
    val denseCentered = centered.toDense(StoragePolicy.AllowDense).toOption.get
    val weights = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 0.0),
        Vector(0.0, 1.0),
        Vector(1.0, 1.0)
      )
    )

    assertMatrixClose(
      centered.rightMultiply(weights).toOption.get,
      GaleNumerics.multiply(denseCentered, weights).toRows,
      1e-12
    )
    assertMatrixClose(centered.crossProduct.toOption.get, GaleNumerics.crossProduct(denseCentered).toRows, 1e-12)

    val selected = centered.selectColumns(IndexSet.from(Vector(2, 0), IndexAxis.Feature).toOption.get).toOption.get
    assertEquals(selected.storage, StorageKind.LazyAffine)
    assertMatrixClose(
      selected.toDense(StoragePolicy.AllowDense).toOption.get,
      denseCentered.toRows.map(row => Vector(row(2), row(0))),
      1e-12
    )
  }

  test("row selection preserves sparse and lazy-affine storage") {
    val selectedRows = IndexSet.from(Vector(3, 1), IndexAxis.Row).toOption.get
    val selectedSparse = sparseView.selectRows(selectedRows).toOption.get
    val centered = PreprocessSpec.Center.fit(sparseView).toOption.get.transform(sparseView).toOption.get
    val selectedCentered = centered.selectRows(selectedRows).toOption.get

    assertEquals(selectedSparse.storage, StorageKind.Sparse)
    assertEquals(selectedCentered.storage, StorageKind.LazyAffine)
    assertMatrixClose(selectedSparse.toDense(StoragePolicy.AllowDense).toOption.get, Vector(denseRows(3), denseRows(1)), 1e-12)
    val transformed = denseTransform(Vector(1.0, 1.0, 1.0), Vector(-1.25, -2.25, -1.75))
    assertMatrixClose(
      selectedCentered.toDense(StoragePolicy.AllowDense).toOption.get,
      Vector(transformed(3), transformed(1)),
      1e-12
    )
  }

  test("MatrixView vector kernels match matrix multiplication without materializing a column matrix") {
    val centered = PreprocessSpec.Center.fit(sparseView).toOption.get.transform(sparseView).toOption.get
    val input = DVec.fromSeq(Vector(2.0, -1.0, 0.5))
    val output = MutableDVec.zeros(centered.rows)
    val transposeInput = DVec.fromSeq(Vector(1.0, -2.0, 0.5, 3.0))
    val transposeOutput = MutableDVec.zeros(centered.cols)

    assert(centered.multiplyVector(input, output).isRight)
    assert(centered.transposeMultiplyVector(transposeInput, transposeOutput).isRight)
    val expected = centered.rightMultiply(
      GaleNumerics.matrixFromRows(Vector.tabulate(input.length)(index => Vector(input(index))))
    ).toOption.get
    val dense = centered.toDense(StoragePolicy.AllowDense).toOption.get
    val expectedTranspose = dense.t * transposeInput
    var row = 0
    while row < output.length do
      assertEqualsDouble(output(row), expected(row, 0), 1e-12)
      row += 1
    var column = 0
    while column < transposeOutput.length do
      assertEqualsDouble(transposeOutput(column), expectedTranspose(column), 1e-12)
      column += 1
  }

  test("standardize scales genuinely varying tiny-magnitude columns to unit variance") {
    val rows = Vector(
      Vector(1e-10 + 5e-13),
      Vector(1e-10 - 5e-13),
      Vector(1e-10 + 5e-13),
      Vector(1e-10 - 5e-13)
    )
    val view = MatrixView.dense(GaleNumerics.matrixFromRows(rows))
    val fitted = PreprocessSpec.Standardize().fit(view).toOption.get
    val out = fitted.transform(view).toOption.get.toDense(StoragePolicy.AllowDense).toOption.get
    val sds = ColumnStats.fromDense(out).flatMap(_.sampleStandardDeviations).toOption.get

    assertEqualsDouble(sds(0), 1.0, 1e-6)
  }

  test("standardize treats constant columns as degenerate and centers them only") {
    val rows = Vector.fill(4)(Vector(3.5, 1e8))
    val view = MatrixView.dense(GaleNumerics.matrixFromRows(rows))
    val fitted = PreprocessSpec.Standardize().fit(view).toOption.get
    val out = fitted.transform(view).toOption.get.toDense(StoragePolicy.AllowDense).toOption.get

    var row = 0
    while row < out.rows do
      assertEqualsDouble(out(row, 0), 0.0, 1e-12)
      assertEqualsDouble(out(row, 1), 0.0, 1e-12)
      row += 1
  }

  test("standardize keeps a unit scale for numerically constant huge-magnitude columns") {
    val rows = Vector.fill(3)(Vector(1e8 + 0.1))
    val view = MatrixView.dense(GaleNumerics.matrixFromRows(rows))
    val fitted = PreprocessSpec.Standardize().fit(view).toOption.get

    fitted match
      case affine: FittedColumnAffine =>
        assertEqualsDouble(affine.scale(0), 1.0, 0.0)
      case other =>
        fail(s"expected FittedColumnAffine, got $other")
  }

  test("fitted preprocessor supports column-subset transform and inverse transform") {
    val fitted = PreprocessSpec.Center.fitInvertible(sparseView).toOption.get
    val columnSet = IndexSet.from(Vector(2, 0), IndexAxis.Feature).toOption.get
    val subset = sparseView.selectColumns(columnSet).toOption.get

    val transformed = fitted.transform(subset, columns = Some(columnSet)).toOption.get
    val restored = fitted.inverseTransform(transformed, columns = Some(columnSet)).toOption.get

    assertMatrixClose(
      transformed.toDense(StoragePolicy.AllowDense).toOption.get,
      denseTransform(Vector(1.0, 1.0, 1.0), Vector(-1.25, -2.25, -1.75)).map(row => Vector(row(2), row(0))),
      1e-12
    )
    assertMatrixClose(
      restored.toDense(StoragePolicy.AllowDense).toOption.get,
      denseRows.map(row => Vector(row(2), row(0))),
      1e-12
    )
  }

  test("preprocessing fit rejects inputs without columns") {
    val view = MatrixView.dense(DMat.zeros(3, 0))

    PreprocessSpec.Pass.fit(view) match
      case Left(MultivarError.InvalidDimension("preprocessing input columns", 0)) => ()
      case other => fail(s"expected InvalidDimension, got $other")
    assert(PreprocessSpec.Standardize().fit(view).isLeft)
  }

  test("column multiplication rejects wrong-length weights") {
    val spec = PreprocessSpec.multiplyColumns(Vector(1.0, 2.0)).toOption.get

    assert(spec.fit(sparseView).swap.toOption.exists(_.message.contains("length 2 != expected 3")))
    assert(PreprocessSpec.multiplyColumns(Vector(1.0, Double.NaN)).isLeft)
  }

  test("requiring invertibility reports a non-invertible scale weight precisely") {
    val fitted = FittedColumnAffine(
      inputCols = 3,
      scale = DVec.fromSeq(Vector(1.0, 0.0, 2.0)),
      shift = MatrixView.zeros(3)
    )

    fitted.requireInvertible match
      case Left(MultivarError.NonInvertibleValue("affine inverse scale", 1, 0.0)) => ()
      case other => fail(s"expected NonInvertibleValue, got $other")
  }

  test("a zero column weight is rejected when the fit must be invertible, not when it is undone") {
    val spec = PreprocessSpec.multiplyColumns(Vector(2.0, 0.0, 1.0)).toOption.get

    assert(spec.fit(sparseView).isRight, "collapsing a column is a valid transform")
    assert(spec.fitInvertible(sparseView).isLeft, "but it cannot be undone")
  }

  test("an invertible fit round-trips and reports its centering and scaling") {
    val fitted = PreprocessSpec.Standardize().fitInvertible(sparseView).toOption.get
    val transformed = fitted.transform(sparseView).toOption.get
    val restored = fitted.inverseTransform(transformed).toOption.get

    assertMatrixClose(restored.toDense(StoragePolicy.AllowDense).toOption.get, denseRows, 1e-12)

    val summary = summaryOf(fitted)
    val stats = sparseView.columnStats.toOption.get

    assertVectorClose(summary.center, stats.means.toOption.get, 1e-12)
    assertVectorClose(summary.scale, stats.sampleStandardDeviations.toOption.get, 1e-12)
  }

  test("dense column-affine inverse matches the view inverseTransform path") {
    val fitted = PreprocessSpec.Standardize().fitInvertible(sparseView).toOption.get
    val working = fitted.transform(sparseView).toOption.get.toDense(StoragePolicy.AllowDense).toOption.get

    val viaView =
      fitted
        .inverseTransform(MatrixView.dense(working), policy = StoragePolicy.AllowDense)
        .flatMap(_.toDense(StoragePolicy.AllowDense))
        .toOption
        .get
    val viaDense = fitted.inverseTransformDense(working).toOption.get

    assertMatrixClose(viaDense, viaView.toRows, 1e-12)
  }

  test("dense column-affine inverse with column selection matches the view path") {
    val fitted = PreprocessSpec.Center.fitInvertible(sparseView).toOption.get
    val columnSet = IndexSet.from(Vector(2, 0), IndexAxis.Feature).toOption.get
    val subset = sparseView.selectColumns(columnSet).toOption.get
    val working = fitted.transform(subset, columns = Some(columnSet)).toOption.get
      .toDense(StoragePolicy.AllowDense)
      .toOption
      .get

    val viaView =
      fitted
        .inverseTransform(MatrixView.dense(working), columns = Some(columnSet), policy = StoragePolicy.AllowDense)
        .flatMap(_.toDense(StoragePolicy.AllowDense))
        .toOption
        .get
    val viaDense = fitted.inverseTransformDense(working, columns = Some(columnSet)).toOption.get

    assertMatrixClose(viaDense, viaView.toRows, 1e-12)
  }

  test("column-affine inverseContributionDense cancels fitted shift in one pass") {
    val fitted = PreprocessSpec.Standardize().fitInvertible(sparseView).toOption.get
    val processed = fitted.transform(sparseView).toOption.get.toDense(StoragePolicy.AllowDense).toOption.get
    val zero = DMat.zeros(processed.rows, processed.cols)

    val original = fitted.inverseTransformDense(processed).toOption.get
    val originalZero = fitted.inverseTransformDense(zero).toOption.get
    val viaDifference = MatrixOps.subtract(original, originalZero)
    val viaContribution = fitted.inverseContributionDense(processed).toOption.get
    val fullInverse = fitted.inverseTransformDense(processed).toOption.get

    assertMatrixClose(viaContribution, viaDifference.toRows, 1e-12)
    assert(!matricesClose(viaContribution, fullInverse, 1e-12), "contribution must cancel shift, not match full inverse")
  }

  private def matricesClose(left: DMat, right: DMat, tol: Double): Boolean =
    left.rows == right.rows && left.cols == right.cols && {
      var row = 0
      var ok = true
      while row < left.rows && ok do
        var col = 0
        while col < left.cols && ok do
          ok = math.abs(left(row, col) - right(row, col)) <= tol
          col += 1
        row += 1
      ok
    }

  test("the variance convention changes the fitted scale") {
    val sample = PreprocessSpec.Standardize().fitInvertible(sparseView).toOption.get
    val population =
      PreprocessSpec.Standardize(VarianceConvention.Population).fitInvertible(sparseView).toOption.get

    val sampleScale = summaryOf(sample).scale
    val populationScale = summaryOf(population).scale
    val rows = sparseView.rows
    val ratio = Math.sqrt((rows - 1).toDouble / rows.toDouble)

    var col = 0
    while col < sampleScale.length do
      assertEqualsDouble(populationScale(col), sampleScale(col) * ratio, 1e-12)
      col += 1
  }

  test("population standard deviations are defined for a single row where sample ones are not") {
    val view = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0, 2.0))))

    PreprocessSpec.Standardize().fit(view) match
      case Left(MultivarError.InsufficientRows("sample standard deviations", 2, 1)) => ()
      case other => fail(s"expected InsufficientRows, got $other")
    assert(PreprocessSpec.Standardize(VarianceConvention.Population).fit(view).isRight)
  }

  test("restrict narrows a fitted preprocessor to selected columns") {
    val fitted = PreprocessSpec.Center.fit(sparseView).toOption.get
    val columnSet = IndexSet.from(Vector(2, 0), IndexAxis.Feature).toOption.get
    val restricted = fitted.restrict(columnSet).toOption.get
    val subset = sparseView.selectColumns(columnSet).toOption.get

    assertEquals(restricted.inputCols, 2)
    assertMatrixClose(
      restricted.transform(subset).toOption.get.toDense(StoragePolicy.AllowDense).toOption.get,
      denseTransform(Vector(1.0, 1.0, 1.0), Vector(-1.25, -2.25, -1.75)).map(row => Vector(row(2), row(0))),
      1e-12
    )
    assert(fitted.restrict(IndexSet.from(Vector(3), IndexAxis.Feature).toOption.get).isLeft)
    assert(fitted.restrict(IndexSet.from(Vector(0), IndexAxis.Row).toOption.get).isLeft)
  }
