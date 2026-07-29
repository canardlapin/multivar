package multivar
package family.spectral

import multivar.core.*

import gale.linalg.Matrix

class MetricCoordinatesSuite extends munit.FunSuite:

  private val tol = 1e-12

  test("table whitening matches manual Q^{1/2} X R^{1/2} for diagonal metrics"):
    val x = Matrix(3, 2)(
      1.0, 2.0,
      3.0, 4.0,
      5.0, 6.0
    )
    val rowMetric = MetricSpec.diagonal(gale.linalg.DVec.fromSeq(Seq(1.0, 4.0, 9.0))).toOption.get
    val featureMetric = MetricSpec.diagonal(gale.linalg.DVec.fromSeq(Seq(4.0, 1.0))).toOption.get
    val coords = MetricCoordinates
      .table(x, rowMetric, featureMetric, DenseSolvers.symmetricEigen, tol)
      .toOption
      .get

    assertEquals(coords.form, MetricCoordinates.Form.Table)
    assertEqualsDouble(coords.matrix(0, 0), 2.0, tol) // 1 * 1 * 2
    assertEqualsDouble(coords.matrix(0, 1), 2.0, tol) // 1 * 2 * 1
    assertEqualsDouble(coords.matrix(1, 0), 12.0, tol) // 2 * 3 * 2
    assertEqualsDouble(coords.matrix(2, 1), 18.0, tol) // 3 * 6 * 1

    val decoded = coords.decode(coords.matrix)
    assertEqualsDouble(decoded(0, 0), 1.0, tol)
    assertEqualsDouble(decoded(2, 1), 6.0, tol)

  test("cross coordinates match G^{-1/2} C R^{1/2} and decode"):
    val cross = Matrix(2, 2)(
      2.0, 0.0,
      0.0, 4.0
    )
    val sourceMetric =
      MetricSpec.diagonal(gale.linalg.DVec.fromSeq(Seq(4.0, 4.0))).toOption.get
    val targetMetric =
      MetricSpec.diagonal(gale.linalg.DVec.fromSeq(Seq(4.0, 1.0))).toOption.get
    val sourceRoots =
      MetricSqrt
        .factor(
          sourceMetric,
          DenseSolvers.symmetricEigen,
          tol,
          StoragePolicy.AllowDense,
          "source metric"
        )
        .toOption
        .get
    val targetRoots =
      MetricSqrt
        .factor(
          targetMetric,
          DenseSolvers.symmetricEigen,
          tol,
          StoragePolicy.AllowDense,
          "target metric"
        )
        .toOption
        .get
    val coords = MetricCoordinates
      .cross(cross, sourceRoots, targetRoots)
      .toOption
      .get

    assertEquals(coords.form, MetricCoordinates.Form.Cross)
    assertEqualsDouble(coords.matrix(0, 0), 2.0, tol) // 0.5 * 2 * 2
    assertEqualsDouble(coords.matrix(1, 1), 2.0, tol) // 0.5 * 4 * 1

    val u = Matrix(2, 1)(1.0, 0.0)
    val v = Matrix(2, 1)(0.0, 1.0)
    assertEqualsDouble(coords.sourceWeights(u)(0, 0), 0.5, tol)
    assertEqualsDouble(coords.targetWeights(v)(1, 0), 1.0, tol)

    val decoded = coords.decode(coords.matrix)
    assertEqualsDouble(decoded(0, 0), 2.0, tol)
    assertEqualsDouble(decoded(1, 1), 4.0, tol)

  test("OrthonormalSubspace projects onto column span without materialising EE^T"):
    val generators = Matrix(3, 2)(
      1.0, 0.0,
      0.0, 0.0,
      0.0, 1.0
    )
    val subspace = OrthonormalSubspace.span(generators, DenseSolvers.svd, tol).toOption.get
    assertEquals(subspace.rank, 2)

    val input = Matrix(3, 1)(3.0, 7.0, 5.0)
    val projected = subspace.project(input).toOption.get
    assertEqualsDouble(projected(0, 0), 3.0, tol)
    assertEqualsDouble(projected(1, 0), 0.0, tol)
    assertEqualsDouble(projected(2, 0), 5.0, tol)

    val split = subspace.split(input).toOption.get
    assertEqualsDouble(split.remainder(1, 0), 7.0, tol)
    assertEqualsDouble(split.focus(1, 0), 0.0, tol)
