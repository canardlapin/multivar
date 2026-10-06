package multivar.inference

import gale.linalg.Matrix
import munit.FunSuite

class MeanProfileSignTestSuite extends FunSuite:
  private def ok[E,A](value: Either[E,A]): A = value.fold(e => fail(e.toString),identity)
  private def matrix(rows: Vector[Vector[Double]]) =
    val b = Matrix.newBuilder(rows.size,rows.head.size)
    for r <- rows.indices; c <- rows.head.indices do b(r,c) = rows(r)(c)
    b.result()
  private def exact(rows: Vector[Vector[Double]]) = ok(ok(MeanProfileSignTest.prepare(matrix(rows))).run(MeanProfileSampling.Exact))
  private def direct(rows: Vector[Vector[Double]],signs: Vector[Int]): Double =
    rows.head.indices.map(c => math.pow(rows.indices.map(r => signs(r)*rows(r)(c)).sum / rows.size,2)).sum

  test("complete orbit has identity, antipode, exact tails and direct signed-mean energies") {
    val rows = Vector(Vector(1.0),Vector(2.0),Vector(3.0))
    val r = exact(rows)
    assertEquals(r.requested,8)
    assertEquals(r.draws.map(_.signs).distinct.size,8)
    assertEquals(r.identityOccurrences,1)
    assertEquals(r.draws.head.signs,Vector(1,1,1))
    assertEquals(r.draws.last.signs,Vector(-1,-1,-1))
    assertEqualsDouble(r.observed,4,1e-14)
    assertEquals(r.exceedances,2)
    assertEqualsDouble(r.pValue.value,.25,1e-15)
    assertEqualsDouble(r.minimumP,.25,1e-15)
    r.draws.foreach(d => assertEqualsDouble(d.statistic,direct(rows,d.signs),1e-14))
  }
  test("Monte Carlo replay of a full orbit retains plus-one rule and external provenance") {
    val prepared = ok(MeanProfileSignTest.prepare(matrix(Vector(Vector(1.0),Vector(2.0),Vector(3.0)))))
    val e = ok(prepared.run(MeanProfileSampling.Exact))
    val r = ok(prepared.run(MeanProfileSampling.ReplayMonteCarlo(e.draws.map(_.signs),"declared uniform replacement stream")))
    assertEquals(r.sampling,"monte-carlo")
    assertEqualsDouble(r.pValue.value,3.0/9,1e-15)
    assertEqualsDouble(r.minimumP,1.0/9,1e-15)
    val a = ok(prepared.run(MeanProfileSampling.MonteCarlo(47,42L)))
    val b = ok(prepared.run(MeanProfileSampling.MonteCarlo(47,42L)))
    assertEquals(a,b)
    val replay = ok(prepared.run(MeanProfileSampling.ReplayMonteCarlo(a.draws.map(_.signs),a.provenance)))
    assertEquals(replay.draws,a.draws)
    assertEquals(replay.pValue,a.pValue)
    assert(a.draws.map(_.signs).distinct.size < a.draws.size)
  }
  test("zero profiles retain all group multiplicities and conservative ties") {
    val z = exact(Vector.fill(3)(Vector(0.0,0.0)))
    assertEquals(z.activeProfiles,0)
    assertEquals(z.groupSize,BigInt(8))
    assertEquals(z.dataOrbitSize,BigInt(1))
    assertEquals(z.exceedances,8)
    assertEqualsDouble(z.pValue.value,1,0)
    assertEqualsDouble(z.minimumP,1,0)
    val inactive = exact(Vector(Vector(1.0),Vector(2.0),Vector(0.0)))
    assertEquals(inactive.activeProfiles,2)
    assertEquals(inactive.dataOrbitSize,BigInt(4))
    assertEqualsDouble(inactive.pValue.value,.5,1e-15)
    assertEqualsDouble(inactive.minimumP,.5,1e-15)
  }
  test("orthogonal coordinates, participant order and measurement scale preserve tails") {
    val rows = Vector(Vector(1.0,3.0),Vector(-2.0,1.0),Vector(.2,-1.0),Vector(4.0,2.0))
    val original = exact(rows)
    val rotated = exact(rows.map(r => Vector((r(0)+r(1))/math.sqrt(2),(r(0)-r(1))/math.sqrt(2))))
    original.draws.zip(rotated.draws).foreach((a,b) => assertEqualsDouble(a.statistic,b.statistic,1e-13))
    assertEquals(original.pValue,rotated.pValue)
    assertEquals(original.pValue,exact(rows.reverse).pValue)
    for scale <- Vector(1e-70,1e70,-7.0) do
      val scaled = exact(rows.map(_.map(_*scale)))
      assertEquals(original.pValue,scaled.pValue)
      assertEqualsDouble(scaled.observed/(scale*scale),original.observed,1e-13)
    original.draws.foreach(d => assertEqualsDouble(d.statistic,direct(rows,d.signs),1e-13))
  }
  test("cancellation and invalid replay preserve completed draws and never return a p-value") {
    val prepared = ok(MeanProfileSignTest.prepare(matrix(Vector(Vector(1.0),Vector(2.0),Vector(3.0)))))
    var checks = 0
    val stopped = prepared.run(MeanProfileSampling.Exact,() => { checks += 1; checks == 4 }).swap.toOption.get
    assertEquals(stopped.completed.size,2)
    assertEquals(stopped.failed,0)
    assert(stopped.cancelled)
    val failed = prepared.run(MeanProfileSampling.ReplayMonteCarlo(Vector(Vector(1,1,1),Vector(1,0,1)),"explicit signs")).swap.toOption.get
    assertEquals(failed.completed.size,1)
    assertEquals(failed.failed,1)
    assert(!failed.cancelled)
    assert(prepared.run(MeanProfileSampling.MonteCarlo(0,1L)).isLeft)
    val limited = ok(MeanProfileSignTest.prepare(matrix(Vector.fill(4)(Vector(1.0))),limits = MeanProfileLimits(maxDraws = 8)))
    assert(limited.run(MeanProfileSampling.Exact).isLeft)
    assert(MeanProfileSignTest.prepare(matrix(Vector(Vector(Double.NaN),Vector(1.0)))).isLeft)
    assert(MeanProfileSignTest.prepare(matrix(Vector(Vector(1e300),Vector(1.0)))).isLeft)
  }
  test("declared absolute roundoff allowance is conservative and scales in squared units") {
    val rows = Vector.fill(4)(Vector(1e-16))
    val a = ok(ok(MeanProfileSignTest.prepare(matrix(rows),absoluteTolerance = 1e-30)).run(MeanProfileSampling.Exact))
    val b = ok(ok(MeanProfileSignTest.prepare(matrix(rows.map(_.map(_*1e8))),absoluteTolerance = 1e-14)).run(MeanProfileSampling.Exact))
    assertEqualsDouble(a.pValue.value,1,0)
    assertEquals(a.pValue,b.pValue)
  }
  test("Holm uses the complete declared family and preserves order, ties and endpoints") {
    val raw = Vector(.01,.04,.03).map(p => ok(PValue(p)))
    val holm = ok(Multiplicity.adjust(raw,MultiplicityMethod.Holm)).map(_.value)
    holm.zip(Vector(.03,.06,.06)).foreach((a,b) => assertEqualsDouble(a,b,1e-15))
    assertEquals(ok(Multiplicity.adjust(raw.reverse,MultiplicityMethod.Holm)).map(_.value),holm.reverse)
    val ties = Vector(0.0,.2,.2,1.0).map(p => ok(PValue(p)))
    assertEquals(ok(Multiplicity.adjust(ties,MultiplicityMethod.Holm)).map(_.value),Vector(0.0,.6000000000000001,.6000000000000001,1.0))
    assert(Multiplicity.adjust(Vector.empty,MultiplicityMethod.Holm).isLeft)
  }
