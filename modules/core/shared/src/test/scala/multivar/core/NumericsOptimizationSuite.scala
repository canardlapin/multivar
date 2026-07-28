package multivar
package core

import gale.linalg.DMat
import gale.spectral.EigenOrder
import gale.spectral.EigenSelection

class NumericsOptimizationSuite extends munit.FunSuite:

  test("isApproximatelyIdentity accepts exact and rejects perturbed identity") {
    val identity = DMat.eye(3)
    val perturbed = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 1e-13, 0.0),
        Vector(0.0, 1.0, 0.0),
        Vector(0.0, 0.0, 1.0)
      )
    )
    val rectangular = GaleNumerics.matrixFromRows(Vector(Vector(1.0, 0.0)))

    assert(MatrixOps.isApproximatelyIdentity(identity, 0.0))
    assert(MatrixOps.isApproximatelyIdentity(perturbed, 1e-12))
    assert(!MatrixOps.isApproximatelyIdentity(perturbed, 1e-14))
    assert(!MatrixOps.isApproximatelyIdentity(rectangular, 1e-12))
  }

  test("inverseSquareRoot short-circuits exact identity without a full spectrum") {
    val dimension = 4
    val identity = DMat.eye(dimension)
    val root = MatrixOps.inverseSquareRoot(identity, DenseSolvers.symmetricEigen, tolerance = 0.0).toOption.get

    assertEquals(root.rows, dimension)
    assertEquals(root.cols, dimension)
    var row = 0
    while row < dimension do
      var col = 0
      while col < dimension do
        val expected = if row == col then 1.0 else 0.0
        assertEqualsDouble(root(row, col), expected, 0.0)
        col += 1
      row += 1
  }

  test("counted symmetric eigen returns the leading slice of the full spectrum") {
    val matrix = GaleNumerics.matrixFromRows(
      Vector(
        Vector(5.0, 1.0, 0.0),
        Vector(1.0, 4.0, 0.0),
        Vector(0.0, 0.0, 2.0)
      )
    )
    val full = DenseSolvers.symmetricEigen.decompose(matrix).toOption.get
    val counted =
      DenseSolvers.symmetricEigen
        .decompose(matrix, EigenSelection.Count(2, EigenOrder.LargestAlgebraic))
        .toOption
        .get

    assertEquals(counted.values.length, 2)
    assertEquals(counted.vectors.cols, 2)
    assertEqualsDouble(counted.values(0), full.values(0), 1e-12)
    assertEqualsDouble(counted.values(1), full.values(1), 1e-12)
  }

  test("column-Gram and row-Gram SVD agree on a square rank-one matrix") {
    val input = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(1.0, 2.0),
          Vector(2.0, 4.0)
        )
      )
    )
    assertEquals(input.rows, input.cols)

    val svd = DenseSolvers.svd.decompose(input, ComponentCount(1).toOption.get).toOption.get
    val reconstructed = GaleNumerics.multiply(
      input.rightMultiply(svd.v).toOption.get,
      svd.v.transpose
    )

    assertEqualsDouble(svd.singularValues(0), 5.0, 1e-9)
    assertMatrixClose(reconstructed, input.toDense().toOption.get.toRows, 1e-9)
  }

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
