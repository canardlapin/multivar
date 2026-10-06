package multivar.inference

class StudentizedIntervalSuite extends munit.FunSuite:
  private def ok[A](value: Either[InferenceError, A]): A = value.fold(e => fail(e.message), identity)
  private def estimate(value: Double, se: Double): StudentizedValue = ok(StudentizedValue(value, se))
  private val ids = Vector.tabulate(5)(i => ok(ReplicateId(i + 10)))
  private val values = Vector(7.0, 9.0, 10.0, 11.0, 13.0)
  private def draws: Iterator[StudentizedDraw] = ids.zip(values).iterator.map((id, value) => StudentizedDraw(id, Right(estimate(value, 1))))
  private val options = StudentizedIntervalOptions(.95)

  test("R type-7 intervals on explicit beta/FIR unit/scaled/context/original-score inputs") {
    StudentizedIntervalOracle.cases.foreach { c =>
      val values = c.values.split(",").toVector.map(_.toDouble)
      val ses = c.ses.split(",").toVector.map(_.toDouble)
      val requested = values.indices.map(i => ok(ReplicateId(3*i+7)))
      val draws = requested.indices.iterator.map(i => StudentizedDraw(requested(i), Right(estimate(values(i), ses(i)))))
      val actual = ok(StudentizedInterval.compute(estimate(c.estimate, c.se), requested, draws, options))
      assertEqualsDouble(actual.lower, c.lower, 1e-10, c.name)
      assertEqualsDouble(actual.upper, c.upper, 1e-10, c.name)
    }
    assertEquals(StudentizedIntervalOracle.cases.size, 56)
  }

  test("jackknife SE matches analytic scalar means including duplicate multisets") {
    for x <- Vector(Vector(1.0,2.0,4.0,8.0,11.0), Vector(1.0,1.0,4.0,8.0,8.0), Vector.fill(5)(2.0)) do
      val deletions = x.indices.map(i => x.patch(i, Vector.empty, 1).sum / (x.size - 1))
      val variance = x.map(v => Math.pow(v - x.sum/x.size, 2)).sum / (x.size - 1)
      assertEqualsDouble(ok(JackknifeStandardError.compute(deletions, x.size)), Math.sqrt(variance/x.size), 2e-14)
    assertEqualsDouble(ok(JackknifeStandardError.compute(Vector(1e15, 1e15+1, 1e15+2), 3)), Math.sqrt(4.0/3), 2e-15)
    assertEqualsDouble(ok(JackknifeStandardError.compute(Vector(1e308, -1e308), 2)) / 1e308, 1.0, 3e-15)
    assert(JackknifeStandardError.compute(Vector(1.0,2.0), 3).isLeft)
    assert(JackknifeStandardError.compute(Vector(1.0,Double.NaN), 2).isLeft)
    assert(JackknifeStandardError.compute(Vector(1.0,2.0), 2, () => true).isLeft)
  }

  test("type-7 pivot inversion and translation/positive-scale/sign equivariance") {
    val result = ok(StudentizedInterval.compute(estimate(10, 2), ids, draws, options))
    assertEqualsDouble(result.lower, 4.4, 3e-15)
    assertEqualsDouble(result.upper, 15.6, 3e-15)
    assertEqualsDouble(result.lowerPivotQuantile, -2.8, 1e-15)
    assertEqualsDouble(result.upperPivotQuantile, 2.8, 1e-15)
    assertEquals(result.replicateIds, ids)
    assertEquals(result.quantileConvention, IntervalQuantileConvention.Type7)
    for scale <- Vector(-7.0, 7.0) do
      val shifted = ids.zip(values).iterator.map((id, value) => StudentizedDraw(id, Right(estimate(value*scale+100, Math.abs(scale)))))
      val transformed = ok(StudentizedInterval.compute(estimate(10*scale+100, 2*Math.abs(scale)), ids, shifted, options))
      val expected = Vector(result.lower*scale+100, result.upper*scale+100).sorted
      assertEqualsDouble(transformed.lower, expected.head, 3e-14)
      assertEqualsDouble(transformed.upper, expected.last, 3e-14)
  }

  test("missing, reordered, duplicate, excluded and extra draws cannot yield an interval") {
    val observed = estimate(10, 2)
    assert(StudentizedInterval.compute(observed, ids, draws.take(4), options).isLeft)
    assert(StudentizedInterval.compute(observed, ids, draws.toVector.reverseIterator, options).isLeft)
    assert(StudentizedInterval.compute(observed, ids.updated(4, ids.head), draws, options).isLeft)
    assert(StudentizedInterval.compute(observed, ids, draws ++ Iterator.single(StudentizedDraw(ids.head, Right(observed))), options).isLeft)
    val failure = InferenceError.RankLoss(2, 1)
    val broken = draws.toVector.updated(2, StudentizedDraw(ids(2), Left(failure)))
    assertEquals(StudentizedInterval.compute(observed, ids, broken.iterator, options),
      Left(InferenceError.ReplicateFailure(ids(2), failure)))
  }

  test("zero SE, declared floors, extreme arithmetic, budgets and cancellation fail explicitly") {
    assert(StudentizedValue(1, 0).isLeft)
    assert(StudentizedValue(1, -1).isLeft)
    assert(StudentizedValue(Double.NaN, 1).isLeft)
    assert(StudentizedValue(1, Double.PositiveInfinity).isLeft)
    val observed = estimate(10, 2)
    assert(StudentizedInterval.compute(observed, ids, draws, options.copy(minimumStandardError = 2)).isLeft)
    assert(StudentizedInterval.compute(observed, ids, draws, options.copy(minimumStandardError = 1)).isLeft)
    assert(StudentizedInterval.compute(observed, ids, draws, options.copy(confidenceLevel = 1)).isLeft)
    assert(StudentizedInterval.compute(observed, ids, draws, options.copy(confidenceLevel = Double.MinPositiveValue)).isLeft)
    var consumed = 0
    val lazyDraws = draws.map { d => consumed += 1; d }
    assert(StudentizedInterval.compute(observed, ids, lazyDraws, options.copy(maxWorkingBytes = 1)).isLeft)
    assertEquals(consumed, 0)
    assert(StudentizedInterval.compute(observed, ids, lazyDraws, options, () => true).isLeft)
    assertEquals(consumed, 0)
    val overflowing = ids.iterator.map(id => StudentizedDraw(id, Right(estimate(-1e308, 1e-300))))
    assert(StudentizedInterval.compute(estimate(1e308, 1), ids, overflowing, options).isLeft)
  }
