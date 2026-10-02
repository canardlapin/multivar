package multivar.family.spectral

import gale.linalg.DMat

class SubspaceAgreementSuite extends munit.FunSuite:
  private def right[A](value: Either[?, A]): A = value.fold(error => fail(error.toString), identity)
  private val plane = DMat.dense(3, 2, Vector(1.0, 0.0, 0.0, 1.0, 0.0, 0.0))

  test("a 45 degree within-plane rotation preserves the subspace but not fixed axes"):
    val c = math.sqrt(.5)
    val rotated = DMat.dense(3, 2, Vector(c, -c, c, c, 0.0, 0.0))
    val result = right(SubspaceAgreement.compare(plane, rotated))
    result.principalAngles.foreach(a => assertEqualsDouble(a, 0.0, 3e-8))
    result.fixedAxisAbsoluteCosines.get.foreach(a => assertEqualsDouble(a, c, 1e-14))

  test("independent 3D geometry gives angles zero and pi over three"):
    val other = DMat.dense(3, 2, Vector(1.0, 0.0, 0.0, .5, 0.0, math.sqrt(.75)))
    val result = right(SubspaceAgreement.compare(plane, other))
    assertEqualsDouble(result.principalAngles(0), 0.0, 3e-8)
    assertEqualsDouble(result.principalAngles(1), math.Pi / 3.0, 1e-12)
    assertEqualsDouble(result.fixedAxisAbsoluteCosines.get(1), .5, 1e-14)

  test("signs and column scales preserve coordinates, while a permutation preserves only the subspace"):
    val scaled = DMat.dense(3, 2, Vector(-1e150, 0.0, 0.0, 1e-150, 0.0, 0.0))
    val result = right(SubspaceAgreement.compare(plane, scaled))
    result.principalAngles.foreach(a => assertEqualsDouble(a, 0.0, 3e-8))
    result.fixedAxisAbsoluteCosines.get.foreach(a => assertEqualsDouble(a, 1.0, 1e-14))
    val swapped = right(SubspaceAgreement.compare(plane, DMat.dense(3, 2, Vector(0.0, 1.0, 1.0, 0.0, 0.0, 0.0))))
    swapped.fixedAxisAbsoluteCosines.get.foreach(a => assertEqualsDouble(a, 0.0, 1e-14))

  test("rank changes retain dimensions and cannot masquerade as equal full subspaces"):
    val line = DMat.dense(3, 1, Vector(1.0, 0.0, 0.0))
    val result = right(SubspaceAgreement.compare(plane, line))
    assertEquals((result.leftRank, result.rightRank), (2, 1))
    assertEquals(result.fixedAxisAbsoluteCosines, None)
    assert(SubspaceAgreement.compare(plane, DMat.dense(3, 2, Vector(1.0, 1.0, 0.0, 0.0, 0.0, 0.0))).isLeft)

  test("foreign dimensions, nonfinite values and owned-array budget refuse"):
    assert(SubspaceAgreement.compare(plane, DMat.eye(4)).isLeft)
    assert(SubspaceAgreement.compare(plane, DMat.dense(3, 1, Vector(1.0, Double.NaN, 0.0))).isLeft)
    assert(SubspaceAgreement.compare(plane, plane, maximumElements = 0).isLeft)
