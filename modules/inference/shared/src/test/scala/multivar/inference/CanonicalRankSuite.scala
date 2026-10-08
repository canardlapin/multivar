package multivar.inference

import gale.backend.Backend.given
import gale.linalg.DMat
import resample4s.kernel.{Permutation, Seed}

class CanonicalRankSuite extends munit.FunSuite:
  private def right[A](value: Either[?, A]): A = value.fold(error => fail(error.toString), identity)
  private def walsh(row: Int, mask: Int): Double = if java.lang.Integer.bitCount(row & mask) % 2 == 0 then 1.0 else -1.0
  private val x = DMat.tabulate(16, 3)((i, j) => walsh(i, 1 << j))
  private val y = DMat.tabulate(16, 2)((i, j) => if j == 0 then .8 * walsh(i, 1) + .6 * walsh(i, 8) else .3 * walsh(i, 2) + math.sqrt(.91) * walsh(i, 3))
  private def permutation(rows: Vector[Int]) = right(Permutation.from(IArray.from(rows)))

  test("unequal candidate spaces and stepwise Wilks match independent R QR/SVD fixtures"):
    val problem = right(StepwiseCanonicalRank.from(x, y))
    assertEquals((problem.leftVariables.cols, problem.rightVariables.cols), (3, 2))
    assertEqualsDouble(problem.correlations(0), .8, 1e-12)
    assertEqualsDouble(problem.correlations(1), .3, 1e-12)
    assertEqualsDouble(problem.observedWilks(0), 1.1159619270032231, 1e-12)
    assertEqualsDouble(problem.observedWilks(1), .09431067947124136, 1e-12)
    val values = right(problem.nullStatistics(permutation(Vector.tabulate(16)(i => (i * 5 + 3) % 16))))
    assertEqualsDouble(values(0), 3.4295968561838563, 1e-11)
    assertEqualsDouble(values(1), 2.4079456086518722, 1e-11)

  test("a retained null-space complement prevents the independently demonstrated dropped-complement counterfactual"):
    val problem = right(StepwiseCanonicalRank.from(x, y))
    val rows = Vector(1,5,10,13,9,6,7,4,11,14,2,3,16,15,12,8).map(_ - 1)
    val values = right(problem.nullStatistics(permutation(rows)))
    assertEqualsDouble(values(0), .4892231550030629, 1e-12)
    assertEqualsDouble(values(1), .054936300920663678, 1e-12)
    val invalidLeft = DMat.tabulate(16, 1)((i, _) => problem.leftVariables(rows(i), 1))
    val invalidRight = DMat.tabulate(16, 1)((i, _) => problem.rightVariables(i, 1))
    val invalid = right(StepwiseCanonicalRank.classical(invalidLeft, invalidRight).flatMap(f => StepwiseCanonicalRank.wilks(f.correlations.toVector)))
    assertEqualsDouble(invalid, .027090948245031738, 1e-12)
    assert(values(1) - invalid > .02)

  test("exact rank-zero and rank-one populations retain every candidate hypothesis"):
    for rank <- 0 to 1 do
      val target = DMat.tabulate(16, 2)((i, j) => if j == 0 && rank == 1 then .8 * walsh(i, 1) + .6 * walsh(i, 8) else walsh(i, if j == 0 then 8 else 3))
      val problem = right(StepwiseCanonicalRank.from(x, target))
      assertEquals(problem.candidateRank, 2)
      assertEquals(problem.observedWilks.size, 2)
      assertEqualsDouble(problem.correlations(0), if rank == 1 then .8 else 0.0, 1e-12)
      assertEqualsDouble(problem.correlations(1), 0.0, 1e-12)
      assertEqualsDouble(problem.observedWilks(1), 0.0, 1e-12)
      val nulls = problem.nullStatistics(permutation(Vector(1,5,10,13,9,6,7,4,11,14,2,3,16,15,12,8).map(_ - 1)))
      if rank == 0 then assert(nulls.left.toOption.exists(_.isInstanceOf[InferenceError.UnidentifiedCanonicalTail]))
      else
        assertEquals(right(nulls).size, 2)
        assert(right(nulls)(1) > 0.0)

  test("reverse dimension imbalance preserves the independent R null-space fixture"):
    val problem = right(StepwiseCanonicalRank.from(y, x))
    val rows = Vector(1,5,10,13,9,6,7,4,11,14,2,3,16,15,12,8).map(_ - 1)
    val inverse = Vector.tabulate(16)(i => rows.indexOf(i))
    val statistics = right(problem.nullStatistics(permutation(inverse)))
    assertEquals((problem.leftVariables.cols, problem.rightVariables.cols), (2,3))
    assertEqualsDouble(statistics(0), .4892231550030629, 1e-11)
    assertEqualsDouble(statistics(1), .054936300920663678, 1e-11)

  test("nuisance residual coordinates preserve the independent canonical roots without recentering the residual basis"):
    val nuisance = DMat.tabulate(16, 2)((i, j) => if j == 0 then 1.0 else walsh(i, 12))
    val basis = right(CanonicalResidualBasis.from(nuisance))
    assertEquals(basis.nuisanceRank, 2)
    val problem = right(StepwiseCanonicalRank.from(right(basis.project(x)), right(basis.project(y))))
    assertEqualsDouble(problem.correlations(0), .8, 1e-12)
    assertEqualsDouble(problem.correlations(1), .3, 1e-12)

  test("Theil basis matches an independent inverse-half formula for an intercept-only design"):
    val n = 8
    val selection = Vector.range(0, 7)
    val basis = right(CanonicalResidualBasis.from(DMat.tabulate(n, 1)((_, _) => 1.0), CanonicalResidualMethod.Theil(selection)))
    val gamma = (math.sqrt(8.0) - 1.0) / 7.0
    for i <- 0 until 8; j <- 0 until 7 do
      val expected = (if i == j then 1.0 else 0.0) - .125 + gamma * (if i < 7 then .125 else -.875)
      assertEqualsDouble(basis.matrix(i, j), expected, 1e-11)

  test("invalid residual selections, rank loss, nonfinite inputs and budgets refuse"):
    val nuisance = DMat.tabulate(8, 2)((i, j) => if j == 0 then 1.0 else if i == 7 then 1.0 else 0.0)
    assert(CanonicalResidualBasis.from(nuisance, CanonicalResidualMethod.Theil(Vector(0,1,2,3,4,7))).isLeft)
    assert(CanonicalResidualBasis.from(nuisance, CanonicalResidualMethod.Theil(Vector.fill(6)(0))).isLeft)
    assert(CanonicalResidualBasis.from(nuisance, maximumElements = 0).isLeft)
    val residual = right(CanonicalResidualBasis.from(nuisance))
    val input = DMat.tabulate(8, 2)((i, j) => walsh(i, 1 << j))
    assert(residual.project(input, maximumOutputElements = 11).isLeft)
    assert(residual.project(input, maximumOutputElements = 12).isRight)
    assert(StepwiseCanonicalRank.from(x, y, maximumElements = 0).isLeft)
    assert(StepwiseCanonicalRank.from(DMat.tabulate(16, 2)((i, _) => walsh(i, 1)), y).isLeft)
    assert(StepwiseCanonicalRank.from(x, DMat.tabulate(16, 2)((_, _) => Double.NaN)).isLeft)
    assert(StepwiseCanonicalRank.from(x, x).isLeft)

  test("fixed B counts inclusive identity ties without early stopping or invented unused draws"):
    val problem = right(StepwiseCanonicalRank.from(x, y))
    val rows = right(RowCount(16))
    val singletonGroups = right(RowPartition.from(rows, Vector.tabulate(16)(i => Vector(i))))
    val result = right(FixedCanonicalRank.run(problem, PermutationAction.withinBlocks(singletonGroups), Seed.fromLong(51L), right(MonteCarloDraws(39)), right(Alpha(.05))))
    assertEquals(result.detectableRank, 0)
    assertEquals(result.identityDraws, 39)
    assertEquals(result.completedCompactFits, 78L)
    assertEquals(result.refitMethod, CanonicalRankRefitMethod.ScoreOrthogonalPermutedTailQrCrossSvdV2)
    result.receipts.foreach: receipt =>
      assertEquals(receipt.allocated.value, 39)
      assertEquals(receipt.consumed.value, 39)
      assertEquals(receipt.exceedances, 39)
      assertEqualsDouble(receipt.pValue.value, 1.0, 0.0)
      assertEquals(receipt.stopReason, MonteCarloStopReason.Exhausted)
      assertEquals(receipt.boundary, None)
      assertEquals(receipt.replicateIds.map(_.value), Vector.range(0, 39))
      assert(receipt.provenance.isInstanceOf[ReplicateProvenance.Deterministic])
    assert(FixedCanonicalRank.run(problem, PermutationAction.unrestricted(rows), Seed.fromLong(51L), right(MonteCarloDraws(39)), right(Alpha(.05)), maximumNullValues = 1).isLeft)

  test("closed adjusted p-values are the prefix maximum of complete fixed-B receipts"):
    val problem = right(StepwiseCanonicalRank.from(x, y))
    val result = right(FixedCanonicalRank.run(problem, PermutationAction.unrestricted(right(RowCount(16))), Seed.fromLong(71L), right(MonteCarloDraws(199)), right(Alpha(.05))))
    assertEquals(result.receipts.size, 2)
    assertEqualsDouble(result.adjustedPValues(0).value, result.receipts(0).pValue.value, 0.0)
    assertEqualsDouble(result.adjustedPValues(1).value, math.max(result.receipts(0).pValue.value, result.receipts(1).pValue.value), 0.0)
    result.receipts.foreach: receipt =>
      assertEquals(receipt.consumed.value, 199)
      assertEqualsDouble(receipt.pValue.value, (1.0 + receipt.exceedances) / 200.0, 0.0)
    assertEquals(result.detectableRank, result.adjustedPValues.takeWhile(_.value <= .05).size)

  test("cached step QR bases match fresh full classical fits for every declared permutation"):
    val problem = right(StepwiseCanonicalRank.from(x, y))
    val plans = Vector(
      Vector.range(0, 16),
      Vector.tabulate(16)(i => (i * 5 + 3) % 16),
      Vector(1,5,10,13,9,6,7,4,11,14,2,3,16,15,12,8).map(_ - 1))
    plans.foreach: rows =>
      val reduced = right(problem.nullStatistics(permutation(rows)))
      for step <- 0 until problem.candidateRank do
        val left = DMat.tabulate(16, problem.leftVariables.cols - step)((i, j) => problem.leftVariables(rows(i), j + step))
        val target = DMat.tabulate(16, problem.rightVariables.cols - step)((i, j) => problem.rightVariables(i, j + step))
        val full = right(StepwiseCanonicalRank.classical(left, target).flatMap(f => StepwiseCanonicalRank.wilks(f.correlations.toVector)))
        assertEqualsDouble(reduced(step), full, 1e-11)

  test("initial zero-root tails match fresh full classical fits after non-identity transforms"):
    // The fitted latent-model API deliberately refuses an all-zero spectrum.
    // Identity zero tails are anchored by the analytic zero-Wilks test above;
    // randomized fits here independently exercise its nonzero fit surface.
    val rows = Vector(1,5,10,13,9,6,7,4,11,14,2,3,16,15,12,8).map(_ - 1)
    val plans = Vector(rows, Vector.tabulate(16)(i => rows.indexOf(i)))
    val target = DMat.tabulate(16, 2)((i, j) => if j == 0 then .8 * walsh(i, 1) + .6 * walsh(i, 8) else walsh(i, 3))
    val problem = right(StepwiseCanonicalRank.from(x, target))
    for rows <- plans; step <- 0 until 2 do
      val reduced = right(problem.nullStatistics(permutation(rows)))
      val left = DMat.tabulate(16, problem.leftVariables.cols - step)((i, j) => problem.leftVariables(rows(i), j + step))
      val rightSide = DMat.tabulate(16, problem.rightVariables.cols - step)((i, j) => problem.rightVariables(i, j + step))
      val full = right(StepwiseCanonicalRank.classical(left, rightSide)).correlations.toVector.map(r => -math.log1p(-r*r)).sum
      assertEqualsDouble(reduced(step), full, 1e-11)

  test("distinct non-identity sampling equals the complete tiny group without success conditioning"):
    val left = DMat.dense(3, 1, Vector(1.0,2.0,4.0))
    val target = DMat.dense(3, 1, Vector(3.0,1.0,2.0))
    val problem = right(StepwiseCanonicalRank.from(left,target))
    val action = PermutationAction.unrestricted(right(RowCount(3)))
    val seed = Seed.fromLong(51L)
    val sampling = CanonicalRankSampling.DistinctNonIdentity(1000)
    val result = right(FixedCanonicalRank.run(problem, action, seed, right(MonteCarloDraws(5)), right(Alpha(.05)), sampling = sampling))
    val selected = result.sampledCandidateIds.map(id => right(action.draw(seed,id)).toIArray.toVector)
    assertEquals(selected.distinct.size, 5)
    assert(!selected.contains(Vector(0,1,2)))
    assertEquals(result.candidateDraws, 5 + result.identityDraws + result.duplicateDraws)
    assert(result.candidateDraws > 5)
    val exact = Vector(0,1,2).permutations.map(rows => right(problem.nullStatistics(permutation(rows)))(0)).toVector
    val expected = exact.count(_ >= problem.observedWilks(0)).toDouble / 6.0
    // Independent scalar cross products are [13,11,15,12,17,16]; four of
    // their squares exceed or tie 13^2. No CCA or SVD enters this count.
    assertEqualsDouble(expected, 2.0 / 3.0, 0.0)
    assertEqualsDouble(result.receipts(0).pValue.value, expected, 0.0)
    assertEquals(result.receipts(0).consumed.value, 5)
    assertEquals(result.completedCompactFits, 5L)
    assertEquals(result.sampling, sampling)
    assert(result.observedCanonicalDiagnostics.allConverged)
    assert(result.maximumNullSvdResidual.isFinite)
    assert(result.maximumNullOrthogonalityError.isFinite)
    assert(FixedCanonicalRank.run(problem, action, seed, right(MonteCarloDraws(6)), right(Alpha(.05)), sampling = sampling).isLeft)
    assert(FixedCanonicalRank.run(problem, action, seed, right(MonteCarloDraws(5)), right(Alpha(.05)), sampling = CanonicalRankSampling.DistinctNonIdentity(5)).left.toOption.exists(_.isInstanceOf[InferenceError.CanonicalDrawBudgetExhausted]))
    assert(FixedCanonicalRank.run(problem, action, seed, right(MonteCarloDraws(5)), right(Alpha(.05)), sampling = sampling, maximumTransformElements = 14).isLeft)
