package multivar
package family.spectral

import gale.linalg.DMat
import multivar.core.{GaleNumerics, MultivarError}

class RotationSuite extends munit.FunSuite:
  private def matrix(rows: Vector[Vector[Double]]): DMat =
    GaleNumerics.matrixFromRows(rows)

  private def assertClose(actual: DMat, expected: DMat, tolerance: Double): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assertEqualsDouble(actual(row, col), expected(row, col), tolerance)
        col += 1
      row += 1

  private def identity(size: Int): DMat = DMat.eye(size)

  test("two-row full-rank Kaiser inputs retain a converged receipt") {
    for input <- Vector(matrix(Vector(Vector(1.0, .4), Vector(.1, 1.2))), matrix(Vector(Vector(1.1, .2), Vector(.3, .9)))) do
      val result = FactorRotation.varimax(input).fold(error => fail(error.message), value => value)
      assert(result.receipt.converged, result.receipt.toString)
      assertClose(result.loadings * result.inverseTransform, input, 1e-10)
  }

  test("default Kaiser rank-one varimax and promax stop with finite receipts") {
    val input = matrix(Vector(Vector(1.0), Vector(-2.0), Vector(3.0)))
    val orthogonal = FactorRotation.varimax(input).fold(error => fail(error.message), value => value)
    val oblique = FactorRotation.promax(input).fold(error => fail(error.message), value => value)
    assertClose(orthogonal.loadings * orthogonal.inverseTransform, input, 1e-12)
    assertClose(oblique.loadings * oblique.inverseTransform, input, 1e-12)
    assert(orthogonal.receipt.converged)
    assert(oblique.receipt.baseline.converged)
    assert(orthogonal.receipt.finalCriterion.isFinite)
    assert(orthogonal.receipt.criterionChange.isFinite)
    assert(oblique.receipt.leastSquaresResidual.isFinite)
    assert(oblique.conditionNumber.isFinite)
  }

  test("varimax keeps a rank-one loading direction and records convergence") {
    val input = matrix(Vector(Vector(1.0), Vector(-2.0), Vector(3.0)))
    val result = FactorRotation.varimax(input, RotationOptions(normalization = RotationNormalization.None)).toOption.get

    assertClose(result.loadings, input, 1e-10)
    assertClose(result.transform, identity(1), 1e-12)
    assert(result.receipt.converged)
    assert(result.receipt.finalCriterion >= result.receipt.initialCriterion - 1e-12)
  }

  test("varimax has an orthogonal transform and preserves reconstruction coordinates") {
    val input = matrix(Vector(
      Vector(3.0, 0.2), Vector(2.0, 0.1), Vector(0.1, 2.0), Vector(0.2, 3.0)
    ))
    val result = FactorRotation.varimax(input).toOption.get
    val orthogonality = result.transform.t * result.transform

    assertClose(orthogonality, identity(2), 1e-8)
    assertClose(result.loadings * result.inverseTransform, input, 1e-8)
  }

  test("varimax criterion is non-decreasing from deterministic identity initialization") {
    val input = matrix(Vector(
      Vector(2.0, 0.4), Vector(1.7, -0.2), Vector(-0.3, 1.9), Vector(0.1, 2.3)
    ))
    val result = FactorRotation.varimax(input, RotationOptions(normalization = RotationNormalization.None)).toOption.get

    assert(result.receipt.finalCriterion >= result.receipt.initialCriterion - 1e-10)
  }

  test("varimax criterion is invariant to signed variable permutations") {
    val input = matrix(Vector(
      Vector(2.0, 0.4), Vector(1.7, -0.2), Vector(-0.3, 1.9), Vector(0.1, 2.3)
    ))
    val signedPermutation = matrix(Vector(
      Vector(0.0, 0.0, -1.0, 0.0), Vector(0.0, 1.0, 0.0, 0.0),
      Vector(0.0, 0.0, 0.0, 1.0), Vector(1.0, 0.0, 0.0, 0.0)
    ))
    val first = FactorRotation.varimax(input, RotationOptions(normalization = RotationNormalization.None)).toOption.get
    val second = FactorRotation.varimax(signedPermutation * input, RotationOptions(normalization = RotationNormalization.None)).toOption.get

    assertEqualsDouble(second.receipt.finalCriterion, first.receipt.finalCriterion, 1e-10)
  }

  test("promax reports score correlation from the retained inverse transform") {
    val input = matrix(Vector(
      Vector(3.0, 0.2), Vector(2.0, 0.3), Vector(0.2, 2.0), Vector(0.3, 3.0), Vector(1.5, 1.1)
    ))
    val result = FactorRotation.promax(input).toOption.get
    val covariance = result.inverseTransform * result.inverseTransform.t

    assertClose(result.scoreCorrelation, covariance, 1e-8)
    assertEqualsDouble(covariance(0, 0), 1.0, 1e-8)
    assertEqualsDouble(covariance(1, 1), 1.0, 1e-8)
    assertClose(result.loadings * result.inverseTransform, input, 1e-8)
    assert(result.conditionNumber >= 1.0)
  }

  test("refuses rank deficient and non-finite rotation inputs") {
    val deficient = matrix(Vector(Vector(1.0, 1.0), Vector(2.0, 2.0), Vector(3.0, 3.0)))
    val nonFinite = matrix(Vector(Vector(1.0, Double.NaN), Vector(2.0, 3.0)))

    assert(FactorRotation.promax(deficient).isLeft)
    assert(FactorRotation.varimax(nonFinite).left.exists:
      case MultivarError.NonFiniteValue("varimax input", 1, value) => value.isNaN
      case _ => false
    )
  }

  test("bounded policies report iteration and allocation refusals") {
    val input = matrix(Vector(Vector(2.0, 0.4), Vector(0.3, 2.0), Vector(1.0, 1.2)))

    assert(FactorRotation.varimax(input, RotationOptions(maxIterations = 1, tolerance = 1e-20)).left.exists(_.isInstanceOf[MultivarError.IterationLimitExceeded]))
    assert(FactorRotation.varimax(input, RotationOptions(maxElements = 2)).left.exists(_.isInstanceOf[MultivarError.DimensionOverflow]))
  }

  test("promax refuses a transform above its explicit SVD condition limit") {
    val input = matrix(Vector(
      Vector(4.0, 0.1), Vector(3.0, 0.2), Vector(0.2, 3.0), Vector(0.1, 4.0), Vector(1.0, 1.0)
    ))

    assert(FactorRotation.promax(input, RotationOptions(maxCondition = 1.0)).isLeft)
  }

  test("oblique score normalization uses inverse times inverse transpose") {
    val shear = matrix(Vector(Vector(1.0, 0.5), Vector(0.0, 1.0)))
    val normalized = FactorRotation.normalizeObliqueTransform(shear, 10.0).toOption.get
    val covariance = normalized._2 * normalized._2.t

    assertEqualsDouble(covariance(0, 0), 1.0, 1e-12)
    assertEqualsDouble(covariance(1, 1), 1.0, 1e-12)
    assertEqualsDouble(covariance(0, 1), -1.0 / Math.sqrt(5.0), 1e-12)
    assertClose(normalized._3, covariance, 1e-12)
  }

  test("Kaiser varimax and normalized power-four promax match independent scalar R optimization") {
    val input = matrix(Vector(
      Vector(1.0, 0.4), Vector(0.8, 0.1), Vector(0.9, 0.3),
      Vector(0.3, 0.9), Vector(0.2, 0.8), Vector(0.5, 1.0)
    ))
    val options = RotationOptions(tolerance = 1e-12)
    val varimax = FactorRotation.varimax(input, options).toOption.get
    val promax = FactorRotation.promax(input, options).toOption.get
    val expectedVarimax = matrix(Vector(
      Vector(0.999306322807462677, 0.037240746461734968), Vector(-0.037240746461734968, 0.999306322807462677)
    ))
    val expectedLoadings = matrix(Vector(
      Vector(0.98441002422276869, 0.43696327558472009), Vector(0.79572098359979671, 0.12972322945013426),
      Vector(0.88820346658819593, 0.33330856865780023), Vector(0.26627522502667733, 0.91054791446523686),
      Vector(0.17006866739210458, 0.80689320753831717), Vector(0.46241241494199636, 1.01792669603833019)
    ))
    val expectedPromax = matrix(Vector(
      Vector(1.13824339541532549, -0.30797058407881722), Vector(-0.39241256566527105, 1.16393288820307905)
    ))

    assert(varimax.receipt.converged)
    assertClose(varimax.transform, expectedVarimax, 1e-8)
    assertClose(varimax.loadings, expectedLoadings, 1e-8)
    assertClose(promax.transform, expectedPromax, 1e-8)
    assertEqualsDouble(promax.scoreCorrelation(0, 1), 0.55690905524985512, 1e-8)
  }
